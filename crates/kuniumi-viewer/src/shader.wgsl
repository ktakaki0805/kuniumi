// The globe of kuniumi-viewer (ticket MV-02).
//
// A vertex holds a unit direction and its height above the sea level in
// metres; it is drawn at the radius 1 + e·h/R (planet radii), with
// `relief_scale_per_m` = e/R. The sea pass draws the same mesh at h = 0.
// Colours arrive as sRGB bytes and are decoded to linear light; the render
// target is an sRGB view, so the hardware encodes them again.

struct Uniforms {
    view_proj: mat4x4<f32>,
    // xyz: the direction towards the light, a unit vector.
    light_dir: vec4<f32>,
    // rgb: the sea colour (linear), a: its opacity.
    sea_colour: vec4<f32>,
    relief_scale_per_m: f32,
};

@group(0) @binding(0) var<uniform> u: Uniforms;

struct VertexOut {
    @builtin(position) @invariant clip: vec4<f32>,
    @location(0) world: vec3<f32>,
    @location(1) colour: vec4<f32>,
};

// The light that reaches a face: ambient plus Lambert.
const AMBIENT: f32 = 0.35;

fn srgb_to_linear(c: vec3<f32>) -> vec3<f32> {
    let low = c / 12.92;
    let high = pow((c + vec3<f32>(0.055)) / 1.055, vec3<f32>(2.4));
    return select(high, low, c <= vec3<f32>(0.04045));
}

fn place(direction: vec3<f32>, height_m: f32) -> VertexOut {
    let p = direction * (1.0 + u.relief_scale_per_m * height_m);
    var out: VertexOut;
    out.clip = u.view_proj * vec4<f32>(p, 1.0);
    out.world = p;
    return out;
}

fn shade(colour: vec4<f32>, normal: vec3<f32>) -> vec4<f32> {
    let light = AMBIENT + (1.0 - AMBIENT) * max(dot(normal, u.light_dir.xyz), 0.0);
    return vec4<f32>(colour.rgb * light, colour.a);
}

@vertex
fn vs_terrain(@location(0) position: vec4<f32>, @location(1) colour: vec4<f32>) -> VertexOut {
    var out = place(position.xyz, position.w);
    out.colour = vec4<f32>(srgb_to_linear(colour.rgb), 1.0);
    return out;
}

@vertex
fn vs_sea(@location(0) position: vec4<f32>) -> VertexOut {
    var out = place(position.xyz, 0.0);
    out.colour = u.sea_colour;
    return out;
}

// Terrain: each triangle lit by its own face normal, from the screen-space
// derivatives of the position, turned outwards, so the relief shows.
@fragment
fn fs_terrain(in: VertexOut) -> @location(0) vec4<f32> {
    let n = normalize(cross(dpdx(in.world), dpdy(in.world)));
    let outward = select(-n, n, dot(n, in.world) >= 0.0);
    return shade(in.colour, outward);
}

// Sea: lit as the smooth sphere it stands for.
@fragment
fn fs_sea(in: VertexOut) -> @location(0) vec4<f32> {
    return shade(in.colour, normalize(in.world));
}
