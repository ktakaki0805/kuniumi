//! Seedable, platform-independent pseudo-random numbers.
//!
//! The algorithms here are frozen: changing them changes every generated
//! world, so any change is a world-breaking change (see docs/conventions.md).
//!
//! - [`SplitMix64`] expands a 64-bit seed into generator state.
//! - [`Rng`] is xoshiro256\*\*, the generator used by all simulation stages.
//! - [`Rng::for_stream`] derives an independent generator from the world seed
//!   and a stable stream name such as `"tectonics/rifting"`, so that adding
//!   random draws to one stage never shifts the numbers another stage sees.

/// SplitMix64. Used only to turn a 64-bit seed into xoshiro256\*\* state.
#[derive(Clone, Debug)]
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    /// Creates the generator from a seed.
    pub const fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// Returns the next 64-bit value.
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

/// 64-bit FNV-1a hash. Used to turn stream names into seeds.
pub const fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut hash = 0xCBF2_9CE4_8422_2325_u64;
    let mut i = 0;
    while i < bytes.len() {
        hash ^= bytes[i] as u64;
        hash = hash.wrapping_mul(0x0000_0100_0000_01B3);
        i += 1;
    }
    hash
}

/// The simulation random number generator (xoshiro256\*\*).
#[derive(Clone, Debug)]
pub struct Rng {
    s: [u64; 4],
}

impl Rng {
    /// Creates a generator directly from a 64-bit seed.
    pub fn from_seed(seed: u64) -> Self {
        let mut sm = SplitMix64::new(seed);
        let mut s = [sm.next_u64(), sm.next_u64(), sm.next_u64(), sm.next_u64()];
        // xoshiro must never have an all-zero state.
        if s == [0; 4] {
            s[0] = 1;
        }
        Self { s }
    }

    /// Derives an independent generator for a named stream of a world seed.
    ///
    /// The stream seed is `fnv1a64(seed as little-endian bytes ++ "/" ++ name)`.
    /// Stream names are part of the output contract: renaming a stream is a
    /// world-breaking change.
    pub fn for_stream(seed: u64, name: &str) -> Self {
        let mut hash = 0xCBF2_9CE4_8422_2325_u64;
        for &b in seed.to_le_bytes().iter().chain(b"/").chain(name.as_bytes()) {
            hash ^= b as u64;
            hash = hash.wrapping_mul(0x0000_0100_0000_01B3);
        }
        Self::from_seed(hash)
    }

    /// Returns the next 64-bit value.
    pub fn next_u64(&mut self) -> u64 {
        let s = &mut self.s;
        let result = s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = s[3].rotate_left(45);
        result
    }

    /// Returns a uniformly distributed `f64` in `[0, 1)` with 53 random bits.
    pub fn next_f64(&mut self) -> f64 {
        const SCALE: f64 = 1.0 / (1u64 << 53) as f64;
        (self.next_u64() >> 11) as f64 * SCALE
    }

    /// Returns a uniformly distributed `f64` in `[lo, hi)`.
    pub fn range_f64(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }

    /// Returns a uniformly distributed integer in `[0, n)` without bias.
    ///
    /// # Panics
    ///
    /// Panics if `n == 0`.
    pub fn below(&mut self, n: u64) -> u64 {
        assert!(n > 0, "Rng::below requires n > 0");
        // Lemire's nearly-divisionless method.
        let mut m = u128::from(self.next_u64()) * u128::from(n);
        let mut low = m as u64;
        if low < n {
            let threshold = n.wrapping_neg() % n;
            while low < threshold {
                m = u128::from(self.next_u64()) * u128::from(n);
                low = m as u64;
            }
        }
        (m >> 64) as u64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Golden values were computed with an independent Python implementation.
    // They must never change: these tests run on every CI platform and are the
    // first line of defence for cross-platform reproducibility.

    #[test]
    fn splitmix64_matches_reference() {
        let mut sm = SplitMix64::new(0);
        assert_eq!(sm.next_u64(), 0xE220_A839_7B1D_CDAF);
        assert_eq!(sm.next_u64(), 0x6E78_9E6A_A1B9_65F4);
        assert_eq!(sm.next_u64(), 0x06C4_5D18_8009_454F);
    }

    #[test]
    fn fnv1a64_matches_reference() {
        assert_eq!(fnv1a64(b""), 0xCBF2_9CE4_8422_2325);
        assert_eq!(fnv1a64(b"a"), 0xAF63_DC4C_8601_EC8C);
    }

    #[test]
    fn xoshiro_matches_reference() {
        let mut rng = Rng::from_seed(0);
        assert_eq!(rng.next_u64(), 0x99EC_5F36_CB75_F2B4);
        assert_eq!(rng.next_u64(), 0xBF6E_1F78_4956_452A);
        let mut rng = Rng::from_seed(42);
        assert_eq!(rng.next_u64(), 0x1578_0B2E_0C2E_C716);
        assert_eq!(rng.next_u64(), 0x6104_D986_6D11_3A7E);
    }

    #[test]
    fn streams_match_reference_and_differ() {
        let mut a = Rng::for_stream(42, "tectonics");
        let mut b = Rng::for_stream(42, "climate");
        assert_eq!(a.next_u64(), 0x9EA9_48C9_B609_8820);
        assert_eq!(b.next_u64(), 0xA58F_A9B2_A577_9EE2);
    }

    #[test]
    fn next_f64_matches_reference_bits() {
        let mut rng = Rng::from_seed(42);
        assert_eq!(rng.next_f64().to_bits(), 0x3FB5_780B_2E0C_2EC0);
        assert_eq!(rng.next_f64().to_bits(), 0x3FD8_4136_619B_444E);
    }

    #[test]
    fn below_stays_in_range() {
        let mut rng = Rng::from_seed(7);
        for n in [1_u64, 2, 3, 7, 1000, u64::MAX] {
            for _ in 0..1000 {
                assert!(rng.below(n) < n);
            }
        }
    }
}
