//! PNG encoding of 8-bit RGB and 16-bit grey images, compressed.
//!
//! A thin layer over the `png` crate (D-10) with its default settings:
//! `Balanced` deflate compression and adaptive row filters. The compressed
//! bytes may change with the crate's version; nothing depends on them (the
//! golden tests fingerprint the pixel buffers instead).

use crate::Rgb8;

/// Encodes a 16-bit grayscale image. `pixels` is row-major, top row first.
///
/// # Panics
///
/// Panics if `pixels.len() != width * height` or a dimension is zero.
pub fn encode_gray16(width: u32, height: u32, pixels: &[u16]) -> Vec<u8> {
    check_dims(width, height, pixels.len());
    let data: Vec<u8> = pixels.iter().flat_map(|v| v.to_be_bytes()).collect();
    encode(
        width,
        height,
        ::png::ColorType::Grayscale,
        ::png::BitDepth::Sixteen,
        &data,
    )
}

/// Encodes an 8-bit RGB image. `pixels` is row-major, top row first.
///
/// # Panics
///
/// Panics if `pixels.len() != width * height` or a dimension is zero.
pub fn encode_rgb8(width: u32, height: u32, pixels: &[Rgb8]) -> Vec<u8> {
    check_dims(width, height, pixels.len());
    encode(
        width,
        height,
        ::png::ColorType::Rgb,
        ::png::BitDepth::Eight,
        pixels.as_flattened(),
    )
}

fn check_dims(width: u32, height: u32, len: usize) {
    assert!(width > 0 && height > 0, "PNG dimensions must be non-zero");
    assert_eq!(
        len,
        width as usize * height as usize,
        "pixel count must equal width * height"
    );
}

/// Encodes raw image bytes (rows top first, no filter bytes) into memory.
fn encode(
    width: u32,
    height: u32,
    colour: ::png::ColorType,
    depth: ::png::BitDepth,
    data: &[u8],
) -> Vec<u8> {
    let mut out = Vec::new();
    let mut encoder = ::png::Encoder::new(&mut out, width, height);
    encoder.set_color(colour);
    encoder.set_depth(depth);
    // The dimensions and the data length are checked, and writing into a
    // `Vec` cannot fail, so neither call returns an error.
    let mut writer = encoder
        .write_header()
        .expect("a valid PNG header writes to memory");
    writer
        .write_image_data(data)
        .expect("checked image data writes to memory");
    writer.finish().expect("a complete PNG writes to memory");
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Decodes `bytes` with the `png` crate: width, height, colour type, bit
    /// depth and the raw image bytes.
    fn decode(bytes: &[u8]) -> (u32, u32, ::png::ColorType, ::png::BitDepth, Vec<u8>) {
        let decoder = ::png::Decoder::new(std::io::Cursor::new(bytes));
        let mut reader = decoder.read_info().unwrap();
        let mut buf = vec![0; reader.output_buffer_size().unwrap()];
        let info = reader.next_frame(&mut buf).unwrap();
        buf.truncate(info.buffer_size());
        (
            info.width,
            info.height,
            info.color_type,
            info.bit_depth,
            buf,
        )
    }

    #[test]
    fn rgb8_round_trips() {
        let (w, h) = (37u32, 23u32);
        let pixels: Vec<Rgb8> = (0..w * h)
            .map(|i| [(i * 7) as u8, (i * 13 + 5) as u8, (i / 3) as u8])
            .collect();
        let (dw, dh, colour, depth, data) = decode(&encode_rgb8(w, h, &pixels));
        assert_eq!((dw, dh), (w, h));
        assert_eq!(colour, ::png::ColorType::Rgb);
        assert_eq!(depth, ::png::BitDepth::Eight);
        assert_eq!(data, pixels.as_flattened());
    }

    #[test]
    fn gray16_round_trips() {
        let (w, h) = (5u32, 3u32);
        let pixels: Vec<u16> = (0..w * h).map(|i| (i * 4_681) as u16).collect();
        assert_eq!(pixels[14], 65_534);
        let (dw, dh, colour, depth, data) = decode(&encode_gray16(w, h, &pixels));
        assert_eq!((dw, dh), (w, h));
        assert_eq!(colour, ::png::ColorType::Grayscale);
        assert_eq!(depth, ::png::BitDepth::Sixteen);
        let expected: Vec<u8> = pixels.iter().flat_map(|v| v.to_be_bytes()).collect();
        assert_eq!(data, expected);
    }

    #[test]
    fn smooth_images_are_compressed() {
        let (w, h) = (1024u32, 512u32);
        let pixels: Vec<Rgb8> = (0..h)
            .flat_map(|y| (0..w).map(move |x| [(x / 4) as u8, (y / 2) as u8, 128]))
            .collect();
        let raw = (w * h * 3) as usize;
        let png = encode_rgb8(w, h, &pixels);
        assert!(png.len() * 20 < raw, "{} bytes for {raw} raw", png.len());
    }

    #[test]
    #[should_panic(expected = "width * height")]
    fn rejects_a_wrong_pixel_count() {
        let _ = encode_rgb8(4, 4, &[[0, 0, 0]; 15]);
    }

    #[test]
    #[should_panic(expected = "non-zero")]
    fn rejects_an_empty_image() {
        let _ = encode_gray16(0, 4, &[]);
    }
}
