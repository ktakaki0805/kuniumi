//! Reading the little-endian binary formats of this crate (`.kfield`,
//! `.kgrid`), with errors that name the format and the part of the file.
//!
//! Every array is read in chunks of at most [`CHUNK_BYTES`] and grows as
//! they arrive (a `Vec` at most doubles its capacity), so a count in a
//! corrupt or truncated header never allocates more than about twice what
//! the input holds, plus one chunk.
//!
//! Errors:
//!
//! - [`io::ErrorKind::UnexpectedEof`]: the input ends early, "truncated
//!   `<format>` file: it ends inside the `<part>`";
//! - [`io::ErrorKind::InvalidData`]: the input does not start with the
//!   format's magic number;
//! - [`io::ErrorKind::Unsupported`]: the format version is not the one this
//!   build reads;
//! - any other error of the reader, unchanged.

use std::io::{self, Read};

/// The largest piece read at once (64 KiB).
pub(crate) const CHUNK_BYTES: usize = 1 << 16;

/// A reader of one file of the format `format` (such as `".kfield"`), for
/// the messages.
pub(crate) struct Reader<R> {
    inner: R,
    format: &'static str,
}

impl<R: Read> Reader<R> {
    /// A reader of `inner`, a file of the format `format`.
    pub(crate) fn new(inner: R, format: &'static str) -> Self {
        Self { inner, format }
    }

    /// Reads the magic number and the format version: `magic` and `version`
    /// are the only ones accepted.
    pub(crate) fn header(&mut self, magic: &[u8; 8], version: u32) -> io::Result<()> {
        let mut found = [0u8; 8];
        let got = self.fill(&mut found)?;
        if found[..got] != magic[..got] {
            return Err(self.invalid(format!(
                "not a {} file: it does not start with the magic number {}",
                self.format,
                magic.escape_ascii()
            )));
        }
        if got < found.len() {
            return Err(self.truncated("magic number"));
        }
        let found = self.u32("format version")?;
        if found != version {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                format!(
                    "{} format version {found} is not supported: this build reads only version \
                     {version}",
                    self.format
                ),
            ));
        }
        Ok(())
    }

    /// A little-endian `u32`, the `part` of the file.
    pub(crate) fn u32(&mut self, part: &str) -> io::Result<u32> {
        let mut b = [0u8; 4];
        self.exact(&mut b, part)?;
        Ok(u32::from_le_bytes(b))
    }

    /// A little-endian `u64`, the `part` of the file.
    pub(crate) fn u64(&mut self, part: &str) -> io::Result<u64> {
        let mut b = [0u8; 8];
        self.exact(&mut b, part)?;
        Ok(u64::from_le_bytes(b))
    }

    /// `len` bytes, the `part` of the file.
    pub(crate) fn bytes(&mut self, len: usize, part: &str) -> io::Result<Vec<u8>> {
        let mut out = Vec::with_capacity(len.min(CHUNK_BYTES));
        while out.len() < len {
            let start = out.len();
            let end = start + (len - start).min(CHUNK_BYTES);
            out.resize(end, 0);
            let got = self.fill(&mut out[start..end])?;
            if start + got < end {
                return Err(
                    self.truncated(&format!("{part} (after {} of {len} bytes)", start + got))
                );
            }
        }
        Ok(out)
    }

    /// `count` little-endian `f64` bit patterns, the `part` of the file.
    pub(crate) fn f64s(&mut self, count: usize, part: &str) -> io::Result<Vec<f64>> {
        self.values::<_, 8>(count, part, |b| {
            f64::from_bits(u64::from_le_bytes(b.try_into().expect("8 bytes")))
        })
    }

    /// `count` little-endian `u32`, the `part` of the file.
    pub(crate) fn u32s(&mut self, count: usize, part: &str) -> io::Result<Vec<u32>> {
        self.values::<_, 4>(count, part, |b| {
            u32::from_le_bytes(b.try_into().expect("4 bytes"))
        })
    }

    /// `count` values of `N` bytes each, decoded by `decode`.
    fn values<T, const N: usize>(
        &mut self,
        count: usize,
        part: &str,
        decode: impl Fn(&[u8]) -> T,
    ) -> io::Result<Vec<T>> {
        let per_chunk = CHUNK_BYTES / N;
        let mut out = Vec::with_capacity(count.min(per_chunk));
        let mut buf = vec![0u8; N * count.min(per_chunk)];
        while out.len() < count {
            let k = (count - out.len()).min(per_chunk);
            let got = self.fill(&mut buf[..N * k])?;
            out.extend(buf[..got - got % N].chunks_exact(N).map(&decode));
            if got < N * k {
                return Err(self.truncated(&format!("{part} (after {} of {count})", out.len())));
            }
        }
        Ok(out)
    }

    /// Fills `buf` exactly, the `part` of the file.
    fn exact(&mut self, buf: &mut [u8], part: &str) -> io::Result<()> {
        if self.fill(buf)? < buf.len() {
            return Err(self.truncated(part));
        }
        Ok(())
    }

    /// Reads into `buf` until it is full or the input ends; the number of
    /// bytes read.
    fn fill(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut got = 0;
        while got < buf.len() {
            match self.inner.read(&mut buf[got..]) {
                Ok(0) => break,
                Ok(k) => got += k,
                Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                Err(e) => return Err(e),
            }
        }
        Ok(got)
    }

    /// The error of an input that ends inside `part`.
    fn truncated(&self, part: &str) -> io::Error {
        io::Error::new(
            io::ErrorKind::UnexpectedEof,
            format!("truncated {} file: it ends inside the {part}", self.format),
        )
    }

    /// An [`io::ErrorKind::InvalidData`] error with `message`.
    pub(crate) fn invalid(&self, message: String) -> io::Error {
        io::Error::new(io::ErrorKind::InvalidData, message)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_values_across_chunks() {
        // 20,000 values of 8 bytes span three chunks of 8,192.
        let values: Vec<f64> = (0..20_000).map(|i| f64::from(i) * 0.5 - 7.0).collect();
        let bytes: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
        let mut r = Reader::new(bytes.as_slice(), ".test");
        let back = r.f64s(values.len(), "values").unwrap();
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&back), bits(&values));
        // The same bytes as u32: twice as many.
        let mut r = Reader::new(bytes.as_slice(), ".test");
        let words = r.u32s(2 * values.len(), "words").unwrap();
        assert_eq!(words[0], (-7.0f64).to_bits() as u32);
        assert_eq!(words[1], ((-7.0f64).to_bits() >> 32) as u32);
    }

    #[test]
    fn an_early_end_names_the_part_and_the_count() {
        let bytes = [0u8; 8 * 10_000 + 5];
        let mut r = Reader::new(&bytes[..], ".test");
        let err = r.f64s(20_000, "values").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
        assert_eq!(
            err.to_string(),
            "truncated .test file: it ends inside the values (after 10000 of 20000)"
        );
        let mut r = Reader::new(&bytes[..3], ".test");
        let err = r.u64("count").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
        assert_eq!(
            err.to_string(),
            "truncated .test file: it ends inside the count"
        );
        let mut r = Reader::new(&bytes[..70_000], ".test");
        let err = r.bytes(100_000, "name").unwrap_err();
        assert_eq!(
            err.to_string(),
            "truncated .test file: it ends inside the name (after 70000 of 100000 bytes)"
        );
    }

    #[test]
    fn a_huge_count_in_a_short_input_allocates_little() {
        // The largest count with 16 bytes behind it: the reader fails after
        // one chunk instead of reserving the count.
        let bytes = [0u8; 16];
        let mut r = Reader::new(&bytes[..], ".test");
        let err = r.f64s(usize::MAX, "values").unwrap_err();
        assert_eq!(
            err.to_string(),
            format!(
                "truncated .test file: it ends inside the values (after 2 of {})",
                usize::MAX
            )
        );
        let mut r = Reader::new(&bytes[..], ".test");
        assert!(r.bytes(usize::MAX, "name").is_err());
    }

    #[test]
    fn checks_the_magic_number_and_the_version() {
        let magic = *b"KUNITST\0";
        let mut good = magic.to_vec();
        good.extend(3u32.to_le_bytes());
        Reader::new(good.as_slice(), ".test")
            .header(&magic, 3)
            .unwrap();
        let err = Reader::new(good.as_slice(), ".test")
            .header(&magic, 4)
            .unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::Unsupported);
        assert_eq!(
            err.to_string(),
            ".test format version 3 is not supported: this build reads only version 4"
        );
        // Another magic number, also in a short input.
        for bad in [&b"KUNIFLD\0\x03\0\0\0"[..], b"KX", b"PNG"] {
            let err = Reader::new(bad, ".test").header(&magic, 3).unwrap_err();
            assert_eq!(err.kind(), io::ErrorKind::InvalidData, "{bad:?}");
            assert_eq!(
                err.to_string(),
                "not a .test file: it does not start with the magic number KUNITST\\x00"
            );
        }
        // Every proper prefix of a good header is truncated.
        for len in 0..good.len() {
            let err = Reader::new(&good[..len], ".test")
                .header(&magic, 3)
                .unwrap_err();
            assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof, "{len}");
            let part = if len < 8 {
                "magic number"
            } else {
                "format version"
            };
            assert_eq!(
                err.to_string(),
                format!("truncated .test file: it ends inside the {part}"),
                "{len}"
            );
        }
    }

    /// A reader that returns one byte at a time and is interrupted before
    /// each.
    struct Trickle<'a> {
        bytes: &'a [u8],
        interrupt: bool,
    }

    impl Read for Trickle<'_> {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            self.interrupt = !self.interrupt;
            if self.interrupt {
                return Err(io::ErrorKind::Interrupted.into());
            }
            match (self.bytes.split_first(), buf.first_mut()) {
                (Some((&b, rest)), Some(slot)) => {
                    *slot = b;
                    self.bytes = rest;
                    Ok(1)
                }
                _ => Ok(0),
            }
        }
    }

    #[test]
    fn short_reads_and_interruptions_are_retried() {
        let bytes: Vec<u8> = (0..40u8).collect();
        let mut r = Reader::new(
            Trickle {
                bytes: &bytes,
                interrupt: false,
            },
            ".test",
        );
        assert_eq!(r.u32("a").unwrap(), u32::from_le_bytes([0, 1, 2, 3]));
        assert_eq!(r.u32s(2, "b").unwrap().len(), 2);
        assert_eq!(r.bytes(4, "c").unwrap(), [12, 13, 14, 15]);
        assert_eq!(r.f64s(3, "d").unwrap().len(), 3);
        assert!(r.u32("e").is_err());
    }
}
