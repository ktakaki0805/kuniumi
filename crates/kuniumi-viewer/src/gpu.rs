//! The GPU side: the globe's buffers and pipelines, the panel drawn over
//! it, drawing a frame into a texture view, and screenshots read back from
//! an offscreen texture.
//!
//! Buffers: the vertices (direction and height, `Float32x4`), the vertex
//! colours (`Unorm8x4`, rewritten when the layer changes), the indices
//! (`Uint32`) and one uniform buffer (`shader.wgsl`). Two pipelines share
//! them: the terrain, opaque with depth writes, and the sea, the same mesh
//! at the sea level, translucent, depth-tested without writes. A third
//! pipeline (`overlay.wgsl`, MV-06) draws the panel's quads
//! (`panel::vertices`) last, alpha-blended over the globe, from their own
//! vertex buffer, the font's atlas (`R8Unorm`) and a uniform with the
//! viewport's size. Depth is `Depth32Float`; the colour target is an sRGB
//! view, with 4× multisampling where the adapter supports it.

use core::future::Future;
use core::num::NonZeroU64;
use core::ops::Range;
use core::pin::pin;
use core::task::{Context, Poll, Waker};
use std::sync::mpsc;

use kuniumi_io::Rgb8;
use kuniumi_io::render::BACKGROUND;
use kuniumi_math::det;
use wgpu::util::DeviceExt;

use crate::font;
use crate::mesh::GlobeMesh;
use crate::panel::{self, Quad, VERTEX_BYTES, VERTICES_PER_QUAD};
use crate::view::Frame;

/// The shader of both pipelines.
pub const SHADER: &str = include_str!("shader.wgsl");

/// The shader of the panel's pipeline.
pub const OVERLAY_SHADER: &str = include_str!("overlay.wgsl");

/// The size of the uniform block of `overlay.wgsl` (bytes): one
/// `vec4<f32>`.
const OVERLAY_UNIFORM_SIZE: u64 = 16;

/// The colour of the sea surface (sRGB).
pub const SEA_COLOUR: Rgb8 = [30, 80, 160];

/// The opacity of the sea surface.
pub const SEA_OPACITY: f64 = 0.55;

/// The depth buffer's format.
const DEPTH_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Depth32Float;

/// The format of screenshots: 8-bit RGBA with sRGB encoding.
const SCREENSHOT_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Rgba8UnormSrgb;

/// The size of the uniform block of `shader.wgsl` (bytes): a 4 × 4 matrix,
/// two `vec4<f32>` and one `f32`, padded to a multiple of 16.
const UNIFORM_SIZE: u64 = 112;

/// The multisample count used where supported.
const MSAA: u32 = 4;

/// The size of the panel's first vertex buffer (bytes); it grows to the
/// next power of two when the quads need more.
const OVERLAY_BUFFER_BYTES: u64 = 1 << 16;

/// Runs `future` to completion on this thread. wgpu's native futures (adapter
/// and device requests) are ready when first polled, so this spins at most a
/// few times; no executor crate is needed.
pub fn block_on<F: Future>(future: F) -> F::Output {
    let mut future = pin!(future);
    let mut context = Context::from_waker(Waker::noop());
    loop {
        if let Poll::Ready(output) = future.as_mut().poll(&mut context) {
            return output;
        }
        std::thread::yield_now();
    }
}

/// An sRGB-encoded channel (0–255) as linear light in [0, 1] (IEC 61966-2-1).
fn srgb_to_linear(channel: u8) -> f64 {
    let c = f64::from(channel) / 255.0;
    if c <= 0.04045 {
        c / 12.92
    } else {
        det::pow((c + 0.055) / 1.055, 2.4)
    }
}

/// The device and queue for `adapter`, with the adapter's own limits (the
/// globe's buffers can exceed the default limits on large grids).
///
/// # Errors
///
/// If a buffer of `mesh` exceeds the adapter's largest buffer, or the device
/// cannot be created.
pub fn request_device(
    adapter: &wgpu::Adapter,
    mesh: &GlobeMesh,
) -> Result<(wgpu::Device, wgpu::Queue), String> {
    let limits = adapter.limits();
    let largest = (16 * mesh.vertices.len()).max(4 * mesh.indices.len()) as u64;
    if largest > limits.max_buffer_size {
        return Err(format!(
            "the globe needs a buffer of {largest} bytes, more than this GPU's largest, {} bytes",
            limits.max_buffer_size
        ));
    }
    block_on(adapter.request_device(&wgpu::DeviceDescriptor {
        label: Some("kuniumi-viewer"),
        required_limits: limits,
        ..Default::default()
    }))
    .map_err(|e| format!("cannot open the GPU device: {e}"))
}

/// The multisample count for drawing into `format` with depth.
fn sample_count(adapter: &wgpu::Adapter, format: wgpu::TextureFormat) -> u32 {
    let supports = |f| {
        adapter
            .get_texture_format_features(f)
            .flags
            .sample_count_supported(MSAA)
    };
    if supports(format) && supports(DEPTH_FORMAT) {
        MSAA
    } else {
        1
    }
}

/// The attachments that depend on the viewport size.
struct Targets {
    width_px: u32,
    height_px: u32,
    /// The multisampled colour target, resolved into the frame's view.
    colour: Option<wgpu::TextureView>,
    depth: wgpu::TextureView,
}

/// The globe on the GPU.
pub struct Renderer {
    device: wgpu::Device,
    queue: wgpu::Queue,
    format: wgpu::TextureFormat,
    sample_count: u32,
    terrain: wgpu::RenderPipeline,
    sea: wgpu::RenderPipeline,
    uniforms: wgpu::Buffer,
    bind_group: wgpu::BindGroup,
    vertices: wgpu::Buffer,
    colours: wgpu::Buffer,
    indices: wgpu::Buffer,
    index_count: u32,
    /// The panel's pipeline, its uniform (the viewport's size), its bind
    /// group (the uniform and the atlas) and its vertices.
    overlay: wgpu::RenderPipeline,
    overlay_uniforms: wgpu::Buffer,
    overlay_bind_group: wgpu::BindGroup,
    overlay_vertices: wgpu::Buffer,
    overlay_vertex_count: u32,
    targets: Option<Targets>,
}

/// A vertex buffer for the panel of `size` bytes.
fn overlay_buffer(device: &wgpu::Device, size: u64) -> wgpu::Buffer {
    device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("panel vertices"),
        size,
        usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    })
}

/// The panel's pipeline drawing into `format` with `sample_count` samples,
/// and its bind group for the atlas made on `queue`.
fn overlay_pipeline(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    format: wgpu::TextureFormat,
    sample_count: u32,
) -> (wgpu::RenderPipeline, wgpu::Buffer, wgpu::BindGroup) {
    let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("panel"),
        source: wgpu::ShaderSource::Wgsl(OVERLAY_SHADER.into()),
    });
    let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("panel"),
        entries: &[
            wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: NonZeroU64::new(OVERLAY_UNIFORM_SIZE),
                },
                count: None,
            },
            wgpu::BindGroupLayoutEntry {
                binding: 1,
                visibility: wgpu::ShaderStages::FRAGMENT,
                ty: wgpu::BindingType::Texture {
                    sample_type: wgpu::TextureSampleType::Float { filterable: false },
                    view_dimension: wgpu::TextureViewDimension::D2,
                    multisampled: false,
                },
                count: None,
            },
        ],
    });
    let [width, height] = font::ATLAS_SIZE;
    let atlas = device.create_texture_with_data(
        queue,
        &wgpu::TextureDescriptor {
            label: Some("font atlas"),
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::R8Unorm,
            usage: wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        },
        wgpu::util::TextureDataOrder::LayerMajor,
        &font::atlas(),
    );
    let uniforms = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("panel uniforms"),
        size: OVERLAY_UNIFORM_SIZE,
        usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("panel"),
        layout: &layout,
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: uniforms.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: wgpu::BindingResource::TextureView(
                    &atlas.create_view(&Default::default()),
                ),
            },
        ],
    });
    let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("panel"),
        bind_group_layouts: &[Some(&layout)],
        immediate_size: 0,
    });
    let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some("panel"),
        layout: Some(&pipeline_layout),
        vertex: wgpu::VertexState {
            module: &shader,
            entry_point: Some("vs_overlay"),
            compilation_options: Default::default(),
            buffers: &[Some(wgpu::VertexBufferLayout {
                array_stride: VERTEX_BYTES as u64,
                step_mode: wgpu::VertexStepMode::Vertex,
                attributes: &wgpu::vertex_attr_array![0 => Float32x2, 1 => Float32x2, 2 => Unorm8x4],
            })],
        },
        primitive: wgpu::PrimitiveState {
            topology: wgpu::PrimitiveTopology::TriangleList,
            cull_mode: None,
            ..Default::default()
        },
        // In the globe's pass, so with its depth buffer, but always drawn
        // and never written to it.
        depth_stencil: Some(wgpu::DepthStencilState {
            format: DEPTH_FORMAT,
            depth_write_enabled: Some(false),
            depth_compare: Some(wgpu::CompareFunction::Always),
            stencil: Default::default(),
            bias: Default::default(),
        }),
        multisample: wgpu::MultisampleState {
            count: sample_count,
            ..Default::default()
        },
        fragment: Some(wgpu::FragmentState {
            module: &shader,
            entry_point: Some("fs_overlay"),
            compilation_options: Default::default(),
            targets: &[Some(wgpu::ColorTargetState {
                format,
                blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                write_mask: wgpu::ColorWrites::ALL,
            })],
        }),
        multiview_mask: None,
        cache: None,
    });
    (pipeline, uniforms, bind_group)
}

impl Renderer {
    /// The globe of `mesh` with the vertex colours `colours`, drawn into
    /// views of `format` (an sRGB format).
    pub fn new(
        adapter: &wgpu::Adapter,
        device: wgpu::Device,
        queue: wgpu::Queue,
        format: wgpu::TextureFormat,
        mesh: &GlobeMesh,
        colours: &[[u8; 4]],
    ) -> Self {
        let sample_count = sample_count(adapter, format);
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("globe"),
            source: wgpu::ShaderSource::Wgsl(SHADER.into()),
        });
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("uniforms"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX_FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: NonZeroU64::new(UNIFORM_SIZE),
                },
                count: None,
            }],
        });
        let uniforms = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("uniforms"),
            size: UNIFORM_SIZE,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("uniforms"),
            layout: &bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: uniforms.as_entire_binding(),
            }],
        });
        let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("globe"),
            bind_group_layouts: &[Some(&bind_group_layout)],
            immediate_size: 0,
        });
        let position = wgpu::VertexBufferLayout {
            array_stride: 16,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &wgpu::vertex_attr_array![0 => Float32x4],
        };
        let colour = wgpu::VertexBufferLayout {
            array_stride: 4,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &wgpu::vertex_attr_array![1 => Unorm8x4],
        };
        let pipeline = |label: &str,
                        vertex: &str,
                        fragment: &str,
                        buffers: &[Option<wgpu::VertexBufferLayout>],
                        sea: bool| {
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: Some(label),
                layout: Some(&layout),
                vertex: wgpu::VertexState {
                    module: &shader,
                    entry_point: Some(vertex),
                    compilation_options: Default::default(),
                    buffers,
                },
                primitive: wgpu::PrimitiveState {
                    topology: wgpu::PrimitiveTopology::TriangleList,
                    front_face: wgpu::FrontFace::Ccw,
                    cull_mode: Some(wgpu::Face::Back),
                    ..Default::default()
                },
                depth_stencil: Some(wgpu::DepthStencilState {
                    format: DEPTH_FORMAT,
                    depth_write_enabled: Some(!sea),
                    depth_compare: Some(wgpu::CompareFunction::Less),
                    stencil: Default::default(),
                    bias: Default::default(),
                }),
                multisample: wgpu::MultisampleState {
                    count: sample_count,
                    ..Default::default()
                },
                fragment: Some(wgpu::FragmentState {
                    module: &shader,
                    entry_point: Some(fragment),
                    compilation_options: Default::default(),
                    targets: &[Some(wgpu::ColorTargetState {
                        format,
                        blend: Some(if sea {
                            wgpu::BlendState::ALPHA_BLENDING
                        } else {
                            wgpu::BlendState::REPLACE
                        }),
                        write_mask: wgpu::ColorWrites::ALL,
                    })],
                }),
                multiview_mask: None,
                cache: None,
            })
        };
        let terrain = pipeline(
            "terrain",
            "vs_terrain",
            "fs_terrain",
            &[Some(position.clone()), Some(colour)],
            false,
        );
        let sea = pipeline("sea", "vs_sea", "fs_sea", &[Some(position)], true);
        let vertices = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("vertices"),
            contents: &mesh
                .vertices
                .iter()
                .flatten()
                .flat_map(|x| x.to_ne_bytes())
                .collect::<Vec<u8>>(),
            usage: wgpu::BufferUsages::VERTEX,
        });
        let colours = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("colours"),
            contents: colours.as_flattened(),
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
        });
        let indices = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("indices"),
            contents: &mesh
                .indices
                .iter()
                .flat_map(|k| k.to_ne_bytes())
                .collect::<Vec<u8>>(),
            usage: wgpu::BufferUsages::INDEX,
        });
        let (overlay, overlay_uniforms, overlay_bind_group) =
            overlay_pipeline(&device, &queue, format, sample_count);
        let overlay_vertices = overlay_buffer(&device, OVERLAY_BUFFER_BYTES);
        Self {
            device,
            queue,
            format,
            sample_count,
            terrain,
            sea,
            uniforms,
            bind_group,
            vertices,
            colours,
            indices,
            index_count: mesh.indices.len() as u32,
            overlay,
            overlay_uniforms,
            overlay_bind_group,
            overlay_vertices,
            overlay_vertex_count: 0,
            targets: None,
        }
    }

    /// Replaces the panel's quads (`panel::layout`), drawn over the globe
    /// from the next [`draw`](Self::draw) on; none draws no panel.
    pub fn set_overlay(&mut self, quads: &[Quad]) {
        let bytes = panel::vertices(quads);
        let needed = bytes.len() as u64;
        if needed > self.overlay_vertices.size() {
            self.overlay_vertices = overlay_buffer(&self.device, needed.next_power_of_two());
        }
        if needed > 0 {
            // A multiple of 4 bytes (20 per vertex), as wgpu requires.
            self.queue.write_buffer(&self.overlay_vertices, 0, &bytes);
        }
        self.overlay_vertex_count = (quads.len() * VERTICES_PER_QUAD) as u32;
    }

    /// The device.
    pub fn device(&self) -> &wgpu::Device {
        &self.device
    }

    /// Presents a window's texture after [`draw`](Self::draw).
    pub fn present(&self, texture: wgpu::SurfaceTexture) {
        self.queue.present(texture);
    }

    /// Replaces the vertex colours (a new layer); `colours` has one entry per
    /// vertex, as the colours given to [`new`](Self::new).
    pub fn set_colours(&self, colours: &[[u8; 4]]) {
        self.queue
            .write_buffer(&self.colours, 0, colours.as_flattened());
    }

    /// Paints the vertices `vertices` (one cell's range of
    /// [`GlobeMesh::cell_vertex_offsets`]) in `colour`, opaque: the
    /// highlight of the cell under the cursor, or its colour restored. Each
    /// vertex colour is 4 bytes, so the write is aligned as wgpu requires.
    pub fn set_cell_colour(&self, vertices: Range<u32>, colour: Rgb8) {
        let [r, g, b] = colour;
        let offset = 4 * u64::from(vertices.start);
        let bytes = [r, g, b, 255].repeat(vertices.len());
        self.queue.write_buffer(&self.colours, offset, &bytes);
    }

    /// Draws `frame` into `view` (of the format given to
    /// [`new`](Self::new)), `width_px` × `height_px`, and submits the work.
    pub fn draw(&mut self, view: &wgpu::TextureView, width_px: u32, height_px: u32, frame: &Frame) {
        self.queue
            .write_buffer(&self.uniforms, 0, &uniform_bytes(frame));
        let viewport: Vec<u8> = [width_px as f32, height_px as f32, 0.0, 0.0]
            .iter()
            .flat_map(|x| x.to_ne_bytes())
            .collect();
        self.queue
            .write_buffer(&self.overlay_uniforms, 0, &viewport);
        self.prepare_targets(width_px, height_px);
        let targets = self.targets.as_ref().expect("prepared above");
        let (attachment, resolve_target) = match &targets.colour {
            Some(multisampled) => (multisampled, Some(view)),
            None => (view, None),
        };
        let [r, g, b] = BACKGROUND.map(srgb_to_linear);
        let mut encoder = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("globe"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: attachment,
                    depth_slice: None,
                    resolve_target,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color { r, g, b, a: 1.0 }),
                        store: if resolve_target.is_some() {
                            wgpu::StoreOp::Discard
                        } else {
                            wgpu::StoreOp::Store
                        },
                    },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &targets.depth,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Discard,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
                multiview_mask: None,
            });
            pass.set_bind_group(0, &self.bind_group, &[]);
            pass.set_index_buffer(self.indices.slice(..), wgpu::IndexFormat::Uint32);
            pass.set_vertex_buffer(0, self.vertices.slice(..));
            pass.set_vertex_buffer(1, self.colours.slice(..));
            pass.set_pipeline(&self.terrain);
            pass.draw_indexed(0..self.index_count, 0, 0..1);
            if frame.sea {
                pass.set_pipeline(&self.sea);
                pass.draw_indexed(0..self.index_count, 0, 0..1);
            }
            if self.overlay_vertex_count > 0 {
                pass.set_pipeline(&self.overlay);
                pass.set_bind_group(0, &self.overlay_bind_group, &[]);
                pass.set_vertex_buffer(0, self.overlay_vertices.slice(..));
                pass.draw(0..self.overlay_vertex_count, 0..1);
            }
        }
        self.queue.submit([encoder.finish()]);
    }

    /// Makes the depth and multisampled targets for a viewport anew when its
    /// size changes.
    fn prepare_targets(&mut self, width_px: u32, height_px: u32) {
        let stale = self
            .targets
            .as_ref()
            .is_none_or(|t| (t.width_px, t.height_px) != (width_px, height_px));
        if stale {
            let texture = |label, format| {
                self.device
                    .create_texture(&wgpu::TextureDescriptor {
                        label: Some(label),
                        size: wgpu::Extent3d {
                            width: width_px,
                            height: height_px,
                            depth_or_array_layers: 1,
                        },
                        mip_level_count: 1,
                        sample_count: self.sample_count,
                        dimension: wgpu::TextureDimension::D2,
                        format,
                        usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
                        view_formats: &[],
                    })
                    .create_view(&Default::default())
            };
            self.targets = Some(Targets {
                width_px,
                height_px,
                colour: (self.sample_count > 1)
                    .then(|| texture("multisampled colour", self.format)),
                depth: texture("depth", DEPTH_FORMAT),
            });
        }
    }
}

/// The uniform block of `shader.wgsl` for `frame`, in the machine's byte
/// order.
fn uniform_bytes(frame: &Frame) -> Vec<u8> {
    let [r, g, b] = SEA_COLOUR.map(srgb_to_linear);
    let [lx, ly, lz] = frame.light_dir;
    let mut floats: Vec<f32> = frame.view_proj.as_flattened().to_vec();
    floats.extend([lx, ly, lz, 0.0]);
    floats.extend([r as f32, g as f32, b as f32, SEA_OPACITY as f32]);
    floats.extend([frame.relief_scale_per_m, 0.0, 0.0, 0.0]);
    let bytes: Vec<u8> = floats.iter().flat_map(|x| x.to_ne_bytes()).collect();
    debug_assert_eq!(bytes.len() as u64, UNIFORM_SIZE);
    bytes
}

/// The globe drawn without a window, into a texture: for screenshots and
/// for timing the rendering itself, unpaced by the display.
pub struct Offscreen {
    renderer: Renderer,
    texture: wgpu::Texture,
    view: wgpu::TextureView,
    width_px: u32,
    height_px: u32,
}

impl Offscreen {
    /// The globe of `mesh` with the vertex colours `colours`, drawn into a
    /// texture of `width_px` × `height_px`.
    ///
    /// # Errors
    ///
    /// If no GPU adapter is available or the device cannot be created.
    pub fn new(
        mesh: &GlobeMesh,
        colours: &[[u8; 4]],
        width_px: u32,
        height_px: u32,
    ) -> Result<Self, String> {
        let instance =
            wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle_from_env());
        let adapter = block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            ..Default::default()
        }))
        .map_err(|e| format!("no GPU adapter: {e}"))?;
        let info = adapter.get_info();
        println!(
            "GPU: {} ({:?}, {:?})",
            info.name, info.backend, info.device_type
        );
        let (device, queue) = request_device(&adapter, mesh)?;
        let renderer = Renderer::new(&adapter, device, queue, SCREENSHOT_FORMAT, mesh, colours);
        let texture = renderer.device.create_texture(&wgpu::TextureDescriptor {
            label: Some("offscreen"),
            size: wgpu::Extent3d {
                width: width_px,
                height: height_px,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: SCREENSHOT_FORMAT,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
            view_formats: &[],
        });
        let view = texture.create_view(&Default::default());
        Ok(Self {
            renderer,
            texture,
            view,
            width_px,
            height_px,
        })
    }

    /// The renderer, for changes to the globe's colours and the panel
    /// before a [`draw`](Self::draw).
    pub fn renderer(&mut self) -> &mut Renderer {
        &mut self.renderer
    }

    /// Draws `frame` and waits until the GPU has finished it.
    ///
    /// # Errors
    ///
    /// If the GPU does not finish.
    pub fn draw(&mut self, frame: &Frame) -> Result<(), String> {
        self.renderer
            .draw(&self.view, self.width_px, self.height_px, frame);
        self.wait()
    }

    fn wait(&self) -> Result<(), String> {
        self.renderer
            .device
            .poll(wgpu::PollType::wait_indefinitely())
            .map(|_| ())
            .map_err(|e| format!("the GPU did not finish: {e}"))
    }

    /// The pixels of the last frame drawn, row-major from the top.
    ///
    /// # Errors
    ///
    /// If the texture cannot be read back.
    pub fn pixels(&self) -> Result<Vec<Rgb8>, String> {
        let device = &self.renderer.device;
        // Rows of the copy are padded to wgpu's alignment.
        let row_bytes = 4 * self.width_px;
        let padded = row_bytes.div_ceil(wgpu::COPY_BYTES_PER_ROW_ALIGNMENT)
            * wgpu::COPY_BYTES_PER_ROW_ALIGNMENT;
        let buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("read back"),
            size: u64::from(padded) * u64::from(self.height_px),
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        let mut encoder = device.create_command_encoder(&Default::default());
        encoder.copy_texture_to_buffer(
            self.texture.as_image_copy(),
            wgpu::TexelCopyBufferInfo {
                buffer: &buffer,
                layout: wgpu::TexelCopyBufferLayout {
                    offset: 0,
                    bytes_per_row: Some(padded),
                    rows_per_image: Some(self.height_px),
                },
            },
            self.texture.size(),
        );
        self.renderer.queue.submit([encoder.finish()]);
        let (sender, receiver) = mpsc::channel();
        buffer.map_async(wgpu::MapMode::Read, .., move |result| {
            // The receiver is alive until the result arrives.
            let _ = sender.send(result);
        });
        self.wait()?;
        let failed = |e: &dyn std::fmt::Display| format!("the frame was not read back: {e}");
        receiver
            .recv()
            .map_err(|e| failed(&e))?
            .map_err(|e| failed(&e))?;
        let data = buffer.get_mapped_range(..).map_err(|e| failed(&e))?;
        Ok(data
            .chunks_exact(padded as usize)
            .flat_map(|row| {
                row[..row_bytes as usize]
                    .chunks_exact(4)
                    .map(|p| [p[0], p[1], p[2]])
            })
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use wgpu::naga;

    use super::*;

    #[test]
    fn shader_parses_and_validates() {
        let module = naga::front::wgsl::parse_str(SHADER)
            .unwrap_or_else(|e| panic!("{}", e.emit_to_string(SHADER)));
        naga::valid::Validator::new(
            naga::valid::ValidationFlags::all(),
            naga::valid::Capabilities::empty(),
        )
        .validate(&module)
        .unwrap_or_else(|e| panic!("{e:?}"));
        let mut entry_points: Vec<&str> = module
            .entry_points
            .iter()
            .map(|e| e.name.as_str())
            .collect();
        entry_points.sort_unstable();
        assert_eq!(
            entry_points,
            ["fs_sea", "fs_terrain", "vs_sea", "vs_terrain"]
        );
        // The uniform block has the size the Rust side writes.
        let uniforms = module
            .global_variables
            .iter()
            .find(|(_, v)| v.name.as_deref() == Some("u"))
            .unwrap()
            .1;
        assert_eq!(
            module.types[uniforms.ty].inner.size(module.to_ctx()),
            UNIFORM_SIZE as u32
        );
    }

    #[test]
    fn overlay_shader_parses_and_validates() {
        let module = naga::front::wgsl::parse_str(OVERLAY_SHADER)
            .unwrap_or_else(|e| panic!("{}", e.emit_to_string(OVERLAY_SHADER)));
        naga::valid::Validator::new(
            naga::valid::ValidationFlags::all(),
            naga::valid::Capabilities::empty(),
        )
        .validate(&module)
        .unwrap_or_else(|e| panic!("{e:?}"));
        let mut entry_points: Vec<&str> = module
            .entry_points
            .iter()
            .map(|e| e.name.as_str())
            .collect();
        entry_points.sort_unstable();
        assert_eq!(entry_points, ["fs_overlay", "vs_overlay"]);
        let global = |name: &str| {
            module
                .global_variables
                .iter()
                .find(|(_, v)| v.name.as_deref() == Some(name))
                .unwrap_or_else(|| panic!("no global `{name}`"))
                .1
        };
        // The bindings of `overlay_pipeline`'s bind group layout (review
        // R1-3).
        let binding = |name: &str| global(name).binding;
        assert_eq!(
            binding("viewport"),
            Some(naga::ResourceBinding {
                group: 0,
                binding: 0
            })
        );
        assert_eq!(
            binding("atlas"),
            Some(naga::ResourceBinding {
                group: 0,
                binding: 1
            })
        );
        // The uniform block has the size the Rust side writes.
        let viewport = global("viewport");
        assert_eq!(
            module.types[viewport.ty].inner.size(module.to_ctx()),
            OVERLAY_UNIFORM_SIZE as u32
        );
        // The atlas: a 2D float texture, read without a sampler.
        assert_eq!(
            module.types[global("atlas").ty].inner,
            naga::TypeInner::Image {
                dim: naga::ImageDimension::D2,
                arrayed: false,
                class: naga::ImageClass::Sampled {
                    kind: naga::ScalarKind::Float,
                    multi: false,
                },
            }
        );
        // The vertex inputs match `panel::vertices`: position and texel as
        // two f32 each, the colour as four normalised channels.
        let vs = module
            .entry_points
            .iter()
            .find(|e| e.name == "vs_overlay")
            .unwrap();
        let inputs: Vec<(u32, naga::TypeInner)> = vs
            .function
            .arguments
            .iter()
            .map(|a| match a.binding {
                Some(naga::Binding::Location { location, .. }) => {
                    (location, module.types[a.ty].inner.clone())
                }
                ref other => panic!("{other:?}"),
            })
            .collect();
        let vector = |size| naga::TypeInner::Vector {
            size,
            scalar: naga::Scalar::F32,
        };
        assert_eq!(
            inputs,
            [
                (0, vector(naga::VectorSize::Bi)),
                (1, vector(naga::VectorSize::Bi)),
                (2, vector(naga::VectorSize::Quad)),
            ]
        );
    }

    #[test]
    fn uniform_layout() {
        let frame = Frame {
            view_proj: [
                [1.0, 2.0, 3.0, 4.0],
                [5.0, 6.0, 7.0, 8.0],
                [9.0, 10.0, 11.0, 12.0],
                [13.0, 14.0, 15.0, 16.0],
            ],
            light_dir: [0.0, 0.6, 0.8],
            relief_scale_per_m: 0.25,
            sea: true,
        };
        let bytes = uniform_bytes(&frame);
        let floats: Vec<f32> = bytes
            .chunks_exact(4)
            .map(|b| f32::from_ne_bytes([b[0], b[1], b[2], b[3]]))
            .collect();
        assert_eq!(floats.len(), 28);
        // Columns in order, as WGSL's mat4x4 expects.
        assert_eq!(
            floats[..16],
            (1..=16).map(|k| k as f32).collect::<Vec<_>>()[..]
        );
        assert_eq!(floats[16..20], [0.0, 0.6, 0.8, 0.0]);
        let sea = SEA_COLOUR.map(|c| srgb_to_linear(c) as f32);
        let bits: Vec<u32> = floats[20..23].iter().map(|x| x.to_bits()).collect();
        assert_eq!(bits, sea.map(f32::to_bits));
        assert_eq!(floats[23].to_bits(), (SEA_OPACITY as f32).to_bits());
        assert_eq!(floats[24..], [0.25, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn srgb_decoding() {
        assert_eq!(srgb_to_linear(0).to_bits(), 0.0_f64.to_bits());
        assert_eq!(srgb_to_linear(255).to_bits(), 1.0_f64.to_bits());
        // 10/255 = 0.0392 lies on the linear segment.
        assert_eq!(
            srgb_to_linear(10).to_bits(),
            (10.0_f64 / 255.0 / 12.92).to_bits()
        );
        // The sRGB mid grey 188 is about 50 % linear light (IEC 61966-2-1).
        assert!((srgb_to_linear(188) - 0.5029).abs() < 1e-4);
    }
}
