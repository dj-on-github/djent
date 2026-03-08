//! djent - Longest run CDF computation.
//!
//! Copyright (C) 2017  David Johnston
//! Translated to Zig 0.15.x, 2025
//!
//! Contact. David Johnston dj@deadhat.com

const std = @import("std");

/// Return the probability of the longest run of heads being <= n
/// in a sequence of r uniform coin tosses.
///
/// The original C version used MPFR for arbitrary precision.
/// This Zig version uses f64, sufficient for typical entropy analysis data sizes.
pub fn longestRunCDF(n: u32, r: u64) f64 {
    const nd: f64 = @floatFromInt(n);
    const rd: f64 = @floatFromInt(r);

    const two_to_n_plus_1 = std.math.pow(f64, 2.0, nd + 1.0);

    const top_a = rd + 1.0;
    const bottom_a = two_to_n_plus_1 - nd - 2.0;
    if (bottom_a == 0.0) return 0.0;

    const first = @exp(-top_a / bottom_a);

    const top_b = two_to_n_plus_1 - 1.0;
    const n_plus_2_over_2 = (nd + 2.0) / 2.0;
    const bottom_b = two_to_n_plus_1 - n_plus_2_over_2;
    if (bottom_b == 0.0) return 0.0;

    const second = top_b / bottom_b;

    return first * second;
}
