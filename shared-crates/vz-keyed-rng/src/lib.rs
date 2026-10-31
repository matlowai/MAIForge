//! `vz.keyed-rng/v1`: counter-based keyed hierarchical randomness.
//!
//! Implements ADR 044 section 2 (node GEN-RNG) and the seed path of ADR 042
//! section 10 (amendment A7). There is no stateful generator and no global
//! sequence: a draw is a pure function of (scope path, purpose, counter,
//! attempt), so it never depends on how many tiles exist, which were computed
//! first, or what a neighbour consumed.
//!
//! # Key derivation (normative for v1)
//!
//! - Root key = `blake3::derive_key("vz.keyed-rng/v1", world seed bytes)`. A
//!   `u64` world seed is its 8 little-endian bytes.
//! - Each child key = `blake3::keyed_hash(parent key, segment bytes)`, where the
//!   segment bytes are one tag byte followed by the typed payload: integers
//!   little-endian, strings as a `u32` little-endian byte length then the bytes.
//! - Scope ladder and tags: law set `0x01`, layer (generator digest) `0x02`,
//!   chart `0x03`, tile `0x04` (lod `u8`, ix `i64`, iy `i64`), address `0x05`.
//!   Segments must appear in ladder order; a segment may be skipped (a
//!   chart-wide layer has no tile), and the tags keep skipped paths distinct.
//! - A draw message = purpose `0x06` (string), counter `0x07` (`u64`),
//!   attempt `0x08` (`u32`), in that order. The draw word is the first 8 bytes
//!   (little-endian `u64`) of `blake3::keyed_hash(scope key, message)`. Shaped
//!   draws that need more bits read the following 8-byte words of the same
//!   BLAKE3 XOF output; word 0 is always the plain draw.
//!
//! The attempt sits in the draw message, after purpose and counter (the ADR
//! 042 section 10 ladder order). Default attempt is 0. Changing one
//! candidate's attempt changes only that candidate's draws (ADR 044 G2).
//!
//! # Shaped draws (integer or table methods only)
//!
//! - `uniform_f64`: top 53 bits times 2^-53, in `[0, 1)`; exact.
//! - `below(n)`: Lemire multiply-shift with exact rejection over XOF words.
//! - `range_i64`, `chance`, `pick_weighted`: built on `below`.
//! - `normal_q32`: table method `normal_q32_table_v1` (Q32 fixed point, see
//!   `normal_table.rs` and `tools/gen_normal_table.py`), truncated at about
//!   +/-5.42 sigma. Integer arithmetic only, so bit-stable on every platform.

#![forbid(unsafe_code)]

mod normal_table;
#[doc(hidden)]
pub mod golden;

use std::fmt;

/// The version id recorded in a world's generation pin (`keyed_rng`).
pub const KEYED_RNG_V1: &str = "vz.keyed-rng/v1";
/// The id of the normal-quantile table this version uses.
pub const NORMAL_TABLE_V1: &str = "normal_q32_table_v1";

pub const TAG_LAW_SET: u8 = 0x01;
pub const TAG_LAYER: u8 = 0x02;
pub const TAG_CHART: u8 = 0x03;
pub const TAG_TILE: u8 = 0x04;
pub const TAG_ADDRESS: u8 = 0x05;
pub const TAG_PURPOSE: u8 = 0x06;
pub const TAG_COUNTER: u8 = 0x07;
pub const TAG_ATTEMPT: u8 = 0x08;
pub const TAG_LEGACY_DERIVATION: u8 = 0x40;
pub const TAG_LEGACY_SALT: u8 = 0x41;

/// Rank used after a legacy `[derivation, salt]` path: only draws may follow.
const RANK_TERMINAL: u8 = TAG_ADDRESS;

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct KeyedRngError(pub String);
impl fmt::Display for KeyedRngError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { f.write_str(&self.0) }
}
impl std::error::Error for KeyedRngError {}
type Result<T> = std::result::Result<T, KeyedRngError>;
fn err<T>(s: impl Into<String>) -> Result<T> { Err(KeyedRngError(s.into())) }

/// Printable ASCII, no space, 1..=max bytes.
fn is_graphic(s: &str, max: usize) -> bool {
    !s.is_empty() && s.len() <= max && s.bytes().all(|b| (0x21..=0x7e).contains(&b))
}
/// Purposes: `[a-z0-9._/-]`, 1..=64 bytes.
fn is_purpose(s: &str) -> bool {
    !s.is_empty() && s.len() <= 64
        && s.bytes().all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || matches!(b, b'.' | b'_' | b'-' | b'/'))
}

fn push_str(out: &mut Vec<u8>, s: &str) {
    out.extend_from_slice(&(s.len() as u32).to_le_bytes());
    out.extend_from_slice(s.as_bytes());
}

/// Canonical bytes of one scope segment (tag + typed payload).
pub fn segment_bytes_str(tag: u8, s: &str) -> Vec<u8> {
    let mut v = vec![tag];
    push_str(&mut v, s);
    v
}
pub fn segment_bytes_chart(chart: u32) -> Vec<u8> {
    let mut v = vec![TAG_CHART];
    v.extend_from_slice(&chart.to_le_bytes());
    v
}
pub fn segment_bytes_tile(lod: u8, ix: i64, iy: i64) -> Vec<u8> {
    let mut v = vec![TAG_TILE, lod];
    v.extend_from_slice(&ix.to_le_bytes());
    v.extend_from_slice(&iy.to_le_bytes());
    v
}
/// Canonical draw message: purpose, counter, attempt.
pub fn draw_message(purpose: &str, counter: u64, attempt: u32) -> Vec<u8> {
    let mut v = Vec::with_capacity(1 + 4 + purpose.len() + 1 + 8 + 1 + 4);
    v.push(TAG_PURPOSE);
    push_str(&mut v, purpose);
    v.push(TAG_COUNTER);
    v.extend_from_slice(&counter.to_le_bytes());
    v.push(TAG_ATTEMPT);
    v.extend_from_slice(&attempt.to_le_bytes());
    v
}

/// One node of the scope ladder: a 32-byte key plus the rank of the last
/// segment applied (0 = root). Copyable, immutable, no state.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct Scope {
    key: [u8; 32],
    rank: u8,
}

impl Scope {
    /// Root of a world whose seed is a `u64` (its 8 little-endian bytes).
    pub fn root(seed: u64) -> Scope { Scope::root_bytes(&seed.to_le_bytes()) }

    /// Root over arbitrary seed bytes.
    pub fn root_bytes(seed: &[u8]) -> Scope {
        Scope { key: blake3::derive_key(KEYED_RNG_V1, seed), rank: 0 }
    }

    /// Legacy `SeedPath { seed, derivation, salt }` as the scope
    /// `[derivation, salt]` under its seed (ADR 044 section 2, Compatibility).
    /// Existing artifacts keep their bytes; this mapping is for new code that
    /// must name a legacy path. Only draws may follow it.
    pub fn legacy_seed_path(seed: u64, derivation: &str, salt: &str) -> Scope {
        let root = Scope::root(seed);
        let d = root.child(&segment_bytes_str(TAG_LEGACY_DERIVATION, derivation), RANK_TERMINAL);
        d.child(&segment_bytes_str(TAG_LEGACY_SALT, salt), RANK_TERMINAL)
    }

    fn child(&self, seg: &[u8], rank: u8) -> Scope {
        Scope { key: *blake3::keyed_hash(&self.key, seg).as_bytes(), rank }
    }

    fn next(&self, tag: u8, seg: Vec<u8>, what: &str) -> Result<Scope> {
        if tag <= self.rank {
            return err(format!("{what} segment out of ladder order (after rank {})", self.rank));
        }
        Ok(self.child(&seg, tag))
    }

    /// Law set and ruleset id (one segment), e.g. a law-set digest.
    pub fn law_set(&self, id: &str) -> Result<Scope> {
        if !is_graphic(id, 128) { return err("law set id must be 1..=128 printable ASCII bytes"); }
        self.next(TAG_LAW_SET, segment_bytes_str(TAG_LAW_SET, id), "law set")
    }

    /// Layer = the layer's generator digest (derivation-recipe digest).
    pub fn layer(&self, generator_digest: &str) -> Result<Scope> {
        if !is_graphic(generator_digest, 128) { return err("layer digest must be 1..=128 printable ASCII bytes"); }
        self.next(TAG_LAYER, segment_bytes_str(TAG_LAYER, generator_digest), "layer")
    }

    pub fn chart(&self, chart: u32) -> Result<Scope> {
        self.next(TAG_CHART, segment_bytes_chart(chart), "chart")
    }

    /// Tile key (lod, ix, iy) of TileKeyV2; the chart is its own segment.
    pub fn tile(&self, lod: u8, ix: i64, iy: i64) -> Result<Scope> {
        self.next(TAG_TILE, segment_bytes_tile(lod, ix, iy), "tile")
    }

    /// The layer-defined local address (ADR 044 section 7 `<local>` part).
    pub fn address(&self, local: &str) -> Result<Scope> {
        if !is_graphic(local, 256) { return err("address must be 1..=256 printable ASCII bytes"); }
        self.next(TAG_ADDRESS, segment_bytes_str(TAG_ADDRESS, local), "address")
    }

    pub fn key(&self) -> [u8; 32] { self.key }
    pub fn key_hex(&self) -> String { hex(&self.key) }

    /// One draw. `attempt` is 0 unless a candidate is being retried.
    pub fn draw(&self, purpose: &str, counter: u64, attempt: u32) -> Result<Draw> {
        if !is_purpose(purpose) { return err("purpose must be 1..=64 bytes of [a-z0-9._/-]"); }
        let mut h = blake3::Hasher::new_keyed(&self.key);
        h.update(&draw_message(purpose, counter, attempt));
        Ok(Draw { reader: h.finalize_xof() })
    }
}

pub fn hex(bytes: &[u8]) -> String {
    const H: &[u8; 16] = b"0123456789abcdef";
    let mut s = String::with_capacity(bytes.len() * 2);
    for b in bytes { s.push(H[(b >> 4) as usize] as char); s.push(H[(b & 15) as usize] as char); }
    s
}

/// The XOF output for one (scope, purpose, counter, attempt). Each shaped
/// method consumes the draw, so one draw yields exactly one value.
pub struct Draw {
    reader: blake3::OutputReader,
}

const MASK40: u64 = (1 << 40) - 1;

impl Draw {
    fn word(&mut self) -> u64 {
        let mut b = [0u8; 8];
        self.reader.fill(&mut b);
        u64::from_le_bytes(b)
    }

    /// The plain draw word (first 8 bytes, little-endian).
    pub fn u64(mut self) -> u64 { self.word() }

    /// Uniform in `[0, 1)` from the top 53 bits; exact.
    pub fn uniform_f64(mut self) -> f64 { (self.word() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0) }

    fn below_words(&mut self, n: u64) -> u64 {
        // Lemire (2019) with exact rejection; extra words come from the XOF.
        let mut m = (self.word() as u128) * (n as u128);
        if (m as u64) < n {
            let t = n.wrapping_neg() % n;
            while (m as u64) < t { m = (self.word() as u128) * (n as u128); }
        }
        (m >> 64) as u64
    }

    /// Uniform integer in `[0, n)`, exactly unbiased. `n` must be > 0.
    pub fn below(mut self, n: u64) -> Result<u64> {
        if n == 0 { return err("below(0) is empty"); }
        Ok(self.below_words(n))
    }

    /// Uniform integer in `[lo, hi)`.
    pub fn range_i64(mut self, lo: i64, hi: i64) -> Result<i64> {
        if hi <= lo { return err("range_i64 needs lo < hi"); }
        let span = (hi as i128 - lo as i128) as u64;
        Ok((lo as i128 + self.below_words(span) as i128) as i64)
    }

    /// True with probability exactly `num / den`.
    pub fn chance(mut self, num: u64, den: u64) -> Result<bool> {
        if den == 0 || num > den { return err("chance needs 0 <= num <= den, den > 0"); }
        Ok(self.below_words(den) < num)
    }

    /// Index chosen with probability `weights[i] / sum`; integer only.
    pub fn pick_weighted(mut self, weights: &[u64]) -> Result<usize> {
        let mut total: u64 = 0;
        for w in weights {
            total = total.checked_add(*w).ok_or_else(|| KeyedRngError("weights overflow u64".into()))?;
        }
        if total == 0 { return err("pick_weighted needs a positive total weight"); }
        let mut r = self.below_words(total);
        for (i, w) in weights.iter().enumerate() {
            if r < *w { return Ok(i); }
            r -= *w;
        }
        unreachable!("r < total")
    }

    /// Standard normal in Q32 fixed point (value times 2^32), table method
    /// `normal_q32_table_v1`. Integer only; truncated at about +/-5.42 sigma.
    pub fn normal_q32(mut self) -> i64 {
        use normal_table::{COARSE, TAIL};
        let w0 = self.word();
        let c = (w0 >> 52) as usize;
        let f = (w0 >> 12) & MASK40;
        if c == 0 || c == 4095 {
            let w1 = self.word();
            let s = (w1 >> 52) as usize;
            let f2 = (w1 >> 12) & MASK40;
            let z = TAIL[s] + lerp40(TAIL[s + 1] - TAIL[s], f2);
            if c == 0 { z } else { -z }
        } else {
            COARSE[c - 1] + lerp40(COARSE[c] - COARSE[c - 1], f)
        }
    }

    /// `normal_q32` as f64 (exact conversion: |value| < 2^53).
    pub fn normal_f64(self) -> f64 { self.normal_q32() as f64 / 4_294_967_296.0 }
}

fn lerp40(d: i64, f: u64) -> i64 { ((d as i128 * f as i128) >> 40) as i64 }

/// BLAKE3 over the table's little-endian i64 bytes (COARSE then TAIL); pinned
/// in the golden file so a silently edited table is caught.
pub fn normal_table_digest() -> String {
    let mut h = blake3::Hasher::new();
    for v in normal_table::COARSE.iter().chain(normal_table::TAIL.iter()) { h.update(&v.to_le_bytes()); }
    format!("blake3:{}", h.finalize().to_hex())
}

/// Normal-table sizes, for tests.
pub fn normal_table_lens() -> (usize, usize) { (normal_table::COARSE.len(), normal_table::TAIL.len()) }
