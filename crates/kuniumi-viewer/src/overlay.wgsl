// The panel of kuniumi-viewer (ticket MV-06): quads in pixels over the
// globe, each filled with a glyph of the bitmap font's atlas (`font.rs`).
//
// A vertex holds its position in pixels from the viewport's top-left
// corner, the matching texel coordinate of the atlas and a colour (sRGB
// bytes and alpha). The atlas is read with `textureLoad`, without a sampler:
// quads lie on whole pixels and are an integer number of pixels per texel,
// so each pixel centre falls well inside one texel. The render target is an
// sRGB view, so colours are decoded to linear light here and encoded again
// by the hardware, as in shader.wgsl.

struct Viewport {
    // xy: the viewport's width and height (px).
    size_px: vec4<f32>,
};

@group(0) @binding(0) var<uniform> viewport: Viewport;
@group(0) @binding(1) var atlas: texture_2d<f32>;

struct VertexOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) texel: vec2<f32>,
    @location(1) colour: vec4<f32>,
};

fn srgb_to_linear(c: vec3<f32>) -> vec3<f32> {
    let low = c / 12.92;
    let high = pow((c + vec3<f32>(0.055)) / 1.055, vec3<f32>(2.4));
    return select(high, low, c <= vec3<f32>(0.04045));
}

@vertex
fn vs_overlay(
    @location(0) position_px: vec2<f32>,
    @location(1) texel: vec2<f32>,
    @location(2) colour: vec4<f32>,
) -> VertexOut {
    let size = viewport.size_px.xy;
    var out: VertexOut;
    out.clip = vec4<f32>(
        2.0 * position_px.x / size.x - 1.0,
        1.0 - 2.0 * position_px.y / size.y,
        0.0,
        1.0,
    );
    out.texel = texel;
    out.colour = vec4<f32>(srgb_to_linear(colour.rgb), colour.a);
    return out;
}

// The glyph's coverage (0 or 1) scales the colour's alpha.
@fragment
fn fs_overlay(in: VertexOut) -> @location(0) vec4<f32> {
    let coverage = textureLoad(atlas, vec2<i32>(floor(in.texel)), 0).r;
    return vec4<f32>(in.colour.rgb, in.colour.a * coverage);
}
