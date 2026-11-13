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
//!
//! The file does not store the field's [`FieldKind`]: [`read_field`] returns
//! a [`FieldKind::Quantity`], and a world's `world.json` records which fields
//! are class fields ([`crate::world_dir`]). Neither the encoding nor the
//! [`fingerprint`] depends on the kind.

use std::io::{self, Read, Write};

use kuniumi_math::rng::fnv1a64;
use kuniumi_world::{Field, FieldKind};

use crate::binary::Reader;

/// File magic.
pub const MAGIC: [u8; 8] = *b"KUNIFLD\0";
/// Current format version.
pub const FORMAT_VERSION: u32 = 1;

/// The name of the format in messages.
const FORMAT: &str = ".kfield";

/// Writes `field` in `.kfield` format.
pub fn write_field(mut w: impl Write, field: &Field) -> io::Result<()> {
    w.write_all(&encode(field))
}

/// Reads a field written by [`write_field`], as a [`FieldKind::Quantity`]
/// (the file does not store the kind).
///
/// # Errors
///
/// [`io::ErrorKind::InvalidData`] if the input does not start with
/// [`MAGIC`] or a string is not UTF-8; [`io::ErrorKind::Unsupported`] for a
/// format version other than [`FORMAT_VERSION`];
/// [`io::ErrorKind::UnexpectedEof`] if the input ends early, naming the part
/// of the file. The strings and the values are read in chunks, so lengths
/// that claim more than the input holds allocate no more than about twice
/// the input plus 64 KiB. Bytes after the last value are not read.
pub fn read_field(r: impl Read) -> io::Result<Field> {
    let mut r = Reader::new(r, FORMAT);
    r.header(&MAGIC, FORMAT_VERSION)?;
    let name = read_string(&mut r, "name")?;
    let unit = read_string(&mut r, "unit")?;
    let count = r.u64("value count")?;
    let count = usize::try_from(count).map_err(|_| {
        r.invalid(format!(
            "the {FORMAT} file gives {count} values, more than this platform can address"
        ))
    })?;
    let values = r.f64s(count, "values")?;
    Ok(Field {
        name,
        unit,
        kind: FieldKind::Quantity,
        values,
    })
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

/// A string: its length (`u32`) and its UTF-8 bytes, the `part` of the
/// file.
fn read_string(r: &mut Reader<impl Read>, part: &str) -> io::Result<String> {
    let len = r.u32(&format!("length of the {part}"))?;
    let bytes = r.bytes(len as usize, part)?;
    String::from_utf8(bytes)
        .map_err(|_| r.invalid(format!("the {part} of the {FORMAT} file is not UTF-8")))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> Field {
        Field {
            name: "surface_height".into(),
            unit: "m".into(),
            kind: FieldKind::Quantity,
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
    fn the_kind_is_not_part_of_the_file() {
        let quantity = sample();
        let class = Field {
            kind: FieldKind::Class,
            ..sample()
        };
        assert_eq!(fingerprint(&quantity), fingerprint(&class));
        let mut buf = Vec::new();
        write_field(&mut buf, &class).unwrap();
        assert_eq!(encode(&quantity), buf);
        assert_eq!(
            read_field(buf.as_slice()).unwrap().kind,
            FieldKind::Quantity
        );
    }

    #[test]
    fn rejects_wrong_magic() {
        let err = read_field(&b"NOTAFILE\x01\0\0\0"[..]).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert_eq!(
            err.to_string(),
            "not a .kfield file: it does not start with the magic number KUNIFLD\\x00"
        );
    }

    #[test]
    fn rejects_other_versions() {
        for version in [0u32, 2, u32::MAX] {
            let mut bytes = encode(&sample());
            bytes[8..12].copy_from_slice(&version.to_le_bytes());
            let err = read_field(bytes.as_slice()).unwrap_err();
            assert_eq!(err.kind(), io::ErrorKind::Unsupported);
            assert_eq!(
                err.to_string(),
                format!(
                    ".kfield format version {version} is not supported: this build reads only \
                     version 1"
                )
            );
        }
    }

    #[test]
    fn every_proper_prefix_is_truncated() {
        // 8 magic, 4 version, 4 + 14 name, 4 + 1 unit, 8 count, 6 × 8 values.
        let bytes = encode(&sample());
        assert_eq!(bytes.len(), 91);
        let parts = [
            (0, "magic number"),
            (8, "format version"),
            (12, "length of the name"),
            (16, "name"),
            (30, "length of the unit"),
            (34, "unit"),
            (35, "value count"),
            (43, "values"),
        ];
        for len in 0..bytes.len() {
            let err = read_field(&bytes[..len]).unwrap_err();
            assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof, "{len}: {err}");
            let part = parts
                .iter()
                .rev()
                .find(|(start, _)| *start <= len)
                .unwrap()
                .1;
            let message = err.to_string();
            assert!(
                message.starts_with(&format!(
                    "truncated .kfield file: it ends inside the {part}"
                )),
                "{len}: {message}"
            );
        }
        let err = read_field(&bytes[..43 + 8 * 4 + 3]).unwrap_err();
        assert_eq!(
            err.to_string(),
            "truncated .kfield file: it ends inside the values (after 4 of 6)"
        );
        let err = read_field(&bytes[..20]).unwrap_err();
        assert_eq!(
            err.to_string(),
            "truncated .kfield file: it ends inside the name (after 4 of 14 bytes)"
        );
    }

    #[test]
    fn huge_counts_in_a_short_file_are_truncated() {
        // A name of 4 GiB − 1 bytes and 2⁶⁴ − 1 values, each with a few
        // bytes behind it: the reader fails at the end of the input
        // instead of reserving the claimed size.
        let mut bytes = encode(&sample())[..12].to_vec();
        bytes.extend(u32::MAX.to_le_bytes());
        bytes.extend(b"abc");
        let err = read_field(bytes.as_slice()).unwrap_err();
        assert_eq!(
            err.to_string(),
            "truncated .kfield file: it ends inside the name (after 3 of 4294967295 bytes)"
        );
        let mut bytes = encode(&sample())[..35].to_vec();
        bytes.extend(u64::MAX.to_le_bytes());
        bytes.extend([0; 12]);
        let err = read_field(bytes.as_slice()).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
        assert!(
            err.to_string()
                .starts_with("truncated .kfield file: it ends inside the values (after 1 of "),
            "{err}"
        );
    }

    #[test]
    fn rejects_names_and_units_that_are_not_utf8() {
        for (at, what) in [(16, "name"), (34, "unit")] {
            let mut bytes = encode(&sample());
            bytes[at] = 0xFF;
            let err = read_field(bytes.as_slice()).unwrap_err();
            assert_eq!(err.kind(), io::ErrorKind::InvalidData);
            assert_eq!(
                err.to_string(),
                format!("the {what} of the .kfield file is not UTF-8")
            );
        }
    }
}
