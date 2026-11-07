//! The `.kfield` binary format for one per-cell field.
//!
//! Layout (all integers little-endian):
//!
//! | bytes            | content                                   |
//! | ---------------- | ----------------------------------------- |
//! | 8                | magic `b"KUNIFLD\0"`                      |
//! | 4                | format version (`u32`, currently 1)       |
//! | 4 + n            | field name: length (`u32`) + UTF-8 bytes  |
//! | 4 + n            | SI unit: length (`u32`) + UTF-8 bytes     |
//! | 8                | value count (`u64`)                       |
//! | 8 × count        | values (`f64`, IEEE 754 bit patterns)     |
//!
//! Values round-trip bit-exactly, so [`fingerprint`] of a field is identical
//! on every platform for the same seed and generator version.

use std::io::{self, Read, Write};

use kuniumi_math::rng::fnv1a64;
use kuniumi_world::Field;

/// File magic.
pub const MAGIC: [u8; 8] = *b"KUNIFLD\0";
/// Current format version.
pub const FORMAT_VERSION: u32 = 1;

/// Writes `field` in `.kfield` format.
pub fn write_field(mut w: impl Write, field: &Field) -> io::Result<()> {
    w.write_all(&encode(field))
}

/// Reads a field written by [`write_field`].
pub fn read_field(mut r: impl Read) -> io::Result<Field> {
    let mut magic = [0u8; 8];
    r.read_exact(&mut magic)?;
    if magic != MAGIC {
        return Err(invalid("not a .kfield file"));
    }
    let version = read_u32(&mut r)?;
    if version != FORMAT_VERSION {
        return Err(invalid("unsupported .kfield version"));
    }
    let name = read_string(&mut r)?;
    let unit = read_string(&mut r)?;
    let count = usize::try_from(read_u64(&mut r)?).map_err(|_| invalid("value count too large"))?;
    let mut values = Vec::with_capacity(count.min(1 << 24));
    let mut buf = [0u8; 8];
    for _ in 0..count {
        r.read_exact(&mut buf)?;
        values.push(f64::from_bits(u64::from_le_bytes(buf)));
    }
    Ok(Field { name, unit, values })
}

/// A 64-bit FNV-1a hash of the field's encoded bytes.
///
/// Two fields have the same fingerprint exactly when name, unit and every
/// value's bit pattern are equal (barring hash collisions). Reproducibility
/// tests compare fingerprints against golden values on every CI platform.
pub fn fingerprint(field: &Field) -> u64 {
    fnv1a64(&encode(field))
}

fn encode(field: &Field) -> Vec<u8> {
    let mut out =
        Vec::with_capacity(32 + field.name.len() + field.unit.len() + 8 * field.values.len());
    out.extend_from_slice(&MAGIC);
    out.extend_from_slice(&FORMAT_VERSION.to_le_bytes());
    for s in [&field.name, &field.unit] {
        let len = u32::try_from(s.len()).expect("string longer than 4 GiB");
        out.extend_from_slice(&len.to_le_bytes());
        out.extend_from_slice(s.as_bytes());
    }
    out.extend_from_slice(&(field.values.len() as u64).to_le_bytes());
    for v in &field.values {
        out.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    out
}

fn read_u32(r: &mut impl Read) -> io::Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

fn read_u64(r: &mut impl Read) -> io::Result<u64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(u64::from_le_bytes(b))
}

fn read_string(r: &mut impl Read) -> io::Result<String> {
    let len = read_u32(r)? as usize;
    let mut bytes = vec![0u8; len];
    r.read_exact(&mut bytes)?;
    String::from_utf8(bytes).map_err(|_| invalid("string is not UTF-8"))
}

fn invalid(msg: &'static str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> Field {
        Field {
            name: "surface_height".into(),
            unit: "m".into(),
            values: vec![-3_700.25, 0.0, -0.0, 8_848.0, f64::MIN_POSITIVE, f64::NAN],
        }
    }

    #[test]
    fn round_trips_bit_exactly() {
        let field = sample();
        let mut buf = Vec::new();
        write_field(&mut buf, &field).unwrap();
        let back = read_field(buf.as_slice()).unwrap();
        assert_eq!(back.name, field.name);
        assert_eq!(back.unit, field.unit);
        let bits = |f: &Field| f.values.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&back), bits(&field));
    }

    #[test]
    fn fingerprint_sees_sign_of_zero() {
        let a = sample();
        let mut b = sample();
        b.values[1] = -0.0;
        assert_ne!(fingerprint(&a), fingerprint(&b));
    }

    #[test]
    fn rejects_wrong_magic() {
        let err = read_field(&b"NOTAFILE\x01\0\0\0"[..]).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }
}
