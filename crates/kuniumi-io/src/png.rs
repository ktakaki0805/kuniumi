//! A minimal, dependency-free PNG encoder.
//!
//! Image data is stored with uncompressed ("stored") deflate blocks, so files
//! are larger than usual but the encoder is tiny and obviously correct. Any
//! PNG reader can open the result.

/// 8-bit RGB pixel.
pub type Rgb8 = [u8; 3];

const SIGNATURE: [u8; 8] = [0x89, b'P', b'N', b'G', 0x0D, 0x0A, 0x1A, 0x0A];
const MAX_STORED_BLOCK: usize = 65_535;

/// Encodes a 16-bit grayscale image. `pixels` is row-major, top row first.
///
/// # Panics
///
/// Panics if `pixels.len() != width * height` or a dimension is zero.
pub fn encode_gray16(width: u32, height: u32, pixels: &[u16]) -> Vec<u8> {
    check_dims(width, height, pixels.len());
    let row_bytes = width as usize * 2;
    let mut raw = Vec::with_capacity((row_bytes + 1) * height as usize);
    for row in pixels.chunks_exact(width as usize) {
        raw.push(0); // filter type: None
        for &v in row {
            raw.extend_from_slice(&v.to_be_bytes());
        }
    }
    encode(width, height, 16, 0, &raw)
}

/// Encodes an 8-bit RGB image. `pixels` is row-major, top row first.
///
/// # Panics
///
/// Panics if `pixels.len() != width * height` or a dimension is zero.
pub fn encode_rgb8(width: u32, height: u32, pixels: &[Rgb8]) -> Vec<u8> {
    check_dims(width, height, pixels.len());
    let row_bytes = width as usize * 3;
    let mut raw = Vec::with_capacity((row_bytes + 1) * height as usize);
    for row in pixels.chunks_exact(width as usize) {
        raw.push(0);
        for px in row {
            raw.extend_from_slice(px);
        }
    }
    encode(width, height, 8, 2, &raw)
}

fn check_dims(width: u32, height: u32, len: usize) {
    assert!(width > 0 && height > 0, "PNG dimensions must be non-zero");
    assert_eq!(
        len,
        width as usize * height as usize,
        "pixel count must equal width * height"
    );
}

fn encode(width: u32, height: u32, bit_depth: u8, color_type: u8, raw: &[u8]) -> Vec<u8> {
    let mut ihdr = Vec::with_capacity(13);
    ihdr.extend_from_slice(&width.to_be_bytes());
    ihdr.extend_from_slice(&height.to_be_bytes());
    ihdr.extend_from_slice(&[bit_depth, color_type, 0, 0, 0]);

    let mut out = Vec::new();
    out.extend_from_slice(&SIGNATURE);
    write_chunk(&mut out, b"IHDR", &ihdr);
    write_chunk(&mut out, b"IDAT", &zlib_stored(raw));
    write_chunk(&mut out, b"IEND", &[]);
    out
}

fn write_chunk(out: &mut Vec<u8>, kind: &[u8; 4], data: &[u8]) {
    let len = u32::try_from(data.len()).expect("PNG chunk larger than 4 GiB");
    out.extend_from_slice(&len.to_be_bytes());
    out.extend_from_slice(kind);
    out.extend_from_slice(data);
    let mut crc = Crc32::new();
    crc.update(kind);
    crc.update(data);
    out.extend_from_slice(&crc.finish().to_be_bytes());
}

/// Wraps `data` in a zlib stream made of stored (uncompressed) deflate blocks.
fn zlib_stored(data: &[u8]) -> Vec<u8> {
    let blocks = data.len().div_ceil(MAX_STORED_BLOCK).max(1);
    let mut out = Vec::with_capacity(data.len() + blocks * 5 + 6);
    out.extend_from_slice(&[0x78, 0x01]); // deflate, 32K window, no preset dict
    let mut chunks = data.chunks(MAX_STORED_BLOCK).peekable();
    if chunks.peek().is_none() {
        out.extend_from_slice(&[1, 0, 0, 0xFF, 0xFF]);
    }
    while let Some(chunk) = chunks.next() {
        let is_final = chunks.peek().is_none();
        let len = chunk.len() as u16;
        out.push(u8::from(is_final));
        out.extend_from_slice(&len.to_le_bytes());
        out.extend_from_slice(&(!len).to_le_bytes());
        out.extend_from_slice(chunk);
    }
    out.extend_from_slice(&adler32(data).to_be_bytes());
    out
}

/// CRC-32 (IEEE 802.3), as used by PNG and zlib.
struct Crc32 {
    value: u32,
}

const CRC_TABLE: [u32; 256] = {
    let mut table = [0u32; 256];
    let mut n = 0;
    while n < 256 {
        let mut c = n as u32;
        let mut k = 0;
        while k < 8 {
            c = if c & 1 != 0 {
                0xEDB8_8320 ^ (c >> 1)
            } else {
                c >> 1
            };
            k += 1;
        }
        table[n] = c;
        n += 1;
    }
    table
};

impl Crc32 {
    fn new() -> Self {
        Self { value: 0xFFFF_FFFF }
    }

    fn update(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.value =
                CRC_TABLE[((self.value ^ u32::from(b)) & 0xFF) as usize] ^ (self.value >> 8);
        }
    }

    fn finish(&self) -> u32 {
        self.value ^ 0xFFFF_FFFF
    }
}

fn adler32(data: &[u8]) -> u32 {
    const MOD: u32 = 65_521;
    let (mut a, mut b) = (1u32, 0u32);
    for chunk in data.chunks(5552) {
        for &byte in chunk {
            a += u32::from(byte);
            b += a;
        }
        a %= MOD;
        b %= MOD;
    }
    (b << 16) | a
}

#[cfg(test)]
mod tests {
    use super::*;

    fn crc32(bytes: &[u8]) -> u32 {
        let mut c = Crc32::new();
        c.update(bytes);
        c.finish()
    }

    #[test]
    fn checksums_match_reference_values() {
        assert_eq!(crc32(b"123456789"), 0xCBF4_3926);
        assert_eq!(crc32(b"IEND"), 0xAE42_6082);
        assert_eq!(adler32(b"Wikipedia"), 0x11E6_0398);
    }

    #[test]
    fn png_layout() {
        let png = encode_gray16(2, 1, &[0, 65535]);
        assert_eq!(&png[..8], &SIGNATURE);
        assert_eq!(&png[12..16], b"IHDR");
        assert_eq!(
            &png[png.len() - 12..],
            &[0, 0, 0, 0, b'I', b'E', b'N', b'D', 0xAE, 0x42, 0x60, 0x82]
        );
    }

    #[test]
    fn large_images_span_several_stored_blocks() {
        let pixels = vec![[1, 2, 3]; 300 * 300];
        let png = encode_rgb8(300, 300, &pixels);
        // 300 rows * (1 + 900) bytes = 270_300 bytes of raw data.
        assert!(png.len() > 270_300);
    }
}
