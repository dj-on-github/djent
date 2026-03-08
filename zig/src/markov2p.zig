//! djent / djrandom - Markov 2-parameter model.
//!
//! Copyright (C) 2017  David Johnston
//! Translated to Zig 0.15.x, 2025
//!
//! Contact. David Johnston dj@deadhat.com

const std = @import("std");
const math = std.math;

pub const TransitionPairType = enum {
    equiprobable,
    p000_max,
    p111_max,
    p101_max,
    p010_max,
};

/// Make two probability density functions for all 2^bitwidth symbols.
pub fn makePDF(p01: f64, p10: f64, bitwidth: u32, table0: []f64, table1: []f64) void {
    const size: usize = @as(usize, 1) << @intCast(bitwidth);
    const p00 = 1.0 - p01;
    const p11 = 1.0 - p10;
    var sum0: f64 = 0.0;
    var sum1: f64 = 0.0;

    for (0..size) |x| {
        if (p01 == 0.5 and p10 == 0.5) {
            const uniform = 1.0 / @as(f64, @floatFromInt(size));
            table0[x] = uniform;
            table1[x] = uniform;
            sum0 += uniform;
            sum1 = sum0;
        } else {
            var plist0: f64 = 1.0;
            var plist1: f64 = 1.0;
            if ((x & 0x1) == 0) { plist0 *= p00; plist1 *= p10; } else { plist0 *= p01; plist1 *= p11; }
            for (0..@as(usize, bitwidth - 1)) |i| {
                const bp = (x >> @intCast(i)) & 0x3;
                switch (bp) {
                    0 => { plist0 *= p00; plist1 *= p00; },
                    1 => { plist0 *= p10; plist1 *= p10; },
                    2 => { plist0 *= p01; plist1 *= p01; },
                    3 => { plist0 *= p11; plist1 *= p11; },
                    else => {},
                }
            }
            table0[x] = plist0;
            table1[x] = plist1;
            sum0 += plist0;
            sum1 += plist1;
        }
    }
    const norm_count = @min(size, 256);
    for (0..norm_count) |i| {
        table0[i] /= sum0;
        table1[i] /= sum1;
    }
}

/// Make two cumulative density functions for all 2^bitwidth symbols.
pub fn makeCDF(p01: f64, p10: f64, bitwidth: u32, table0: []f64, table1: []f64) void {
    const size: usize = @as(usize, 1) << @intCast(bitwidth);
    const p00 = 1.0 - p01;
    const p11 = 1.0 - p10;

    for (0..size) |x| {
        if (p01 == 0.5 and p10 == 0.5) {
            const uniform = 1.0 / @as(f64, @floatFromInt(size));
            if (x == 0) { table0[x] = uniform; table1[x] = uniform; } else { table0[x] = table0[x - 1] + uniform; table1[x] = table1[x - 1] + uniform; }
        } else {
            var plist0: f64 = 1.0;
            var plist1: f64 = 1.0;
            if ((x & 0x1) == 0) { plist0 *= p00; plist1 *= p10; } else { plist0 *= p01; plist1 *= p11; }
            for (0..@as(usize, bitwidth - 1)) |i| {
                const bp = (x >> @intCast(i)) & 0x3;
                switch (bp) {
                    0 => { plist0 *= p00; plist1 *= p00; },
                    1 => { plist0 *= p10; plist1 *= p10; },
                    2 => { plist0 *= p01; plist1 *= p01; },
                    3 => { plist0 *= p11; plist1 *= p11; },
                    else => {},
                }
            }
            if (x == 0) { table0[x] = plist0; table1[x] = plist1; } else { table0[x] = table0[x - 1] + plist0; table1[x] = table1[x - 1] + plist1; }
        }
    }
    const max0 = table0[size - 1];
    const max1 = table1[size - 1];
    for (0..size) |i| { table0[i] /= max0; table1[i] /= max1; }
}

/// Compute the symbol probability for the Markov 2-parameter model.
pub fn symbolProb(p01: f64, p10: f64, x: u64, bitwidth: u32) f64 {
    const p00 = 1.0 - p01;
    const p11 = 1.0 - p10;
    const mu = p01 / (p10 + p01);
    const p0 = 1.0 - mu;
    const p1 = mu;

    if (p01 == 0.5 and p10 == 0.5) return 1.0;

    var plist0: f64 = 1.0;
    var plist1: f64 = 1.0;

    if (((x >> @intCast(bitwidth - 1)) & 0x1) == 0) { plist0 *= p00; plist1 *= p10; } else { plist0 *= p01; plist1 *= p11; }

    if (bitwidth >= 2) {
        for (0..@as(usize, bitwidth - 2)) |i| {
            const shift: u6 = @intCast(bitwidth - 2 - @as(u32, @intCast(i)));
            const bp: u2 = @intCast((x >> shift) & 0x3);
            switch (bp) {
                0 => { plist0 *= p00; plist1 *= p00; },
                1 => { plist0 *= p01; plist1 *= p01; },
                2 => { plist0 *= p10; plist1 *= p10; },
                3 => { plist0 *= p11; plist1 *= p11; },
            }
        }
    }
    return (p0 * plist0) + (p1 * plist1);
}

pub fn mostProbableTransitionPair(p01: f64, p10: f64) TransitionPairType {
    const mu = p01 / (p10 + p01);
    const p0 = 1.0 - mu;
    const p1 = mu;
    const p00 = 1.0 - p01;
    const p11 = 1.0 - p10;
    const p010 = p0 * p01 * p10;
    const p101 = p1 * p10 * p01;
    const p000 = p0 * p00 * p00;
    const p111 = p1 * p11 * p11;

    if (p111 >= p000 and p111 >= p101 and p111 >= p010) return .p111_max;
    if (p000 >= p111 and p000 >= p101 and p000 >= p010) return .p000_max;
    if (p101 >= p111 and p101 >= p000 and p101 >= p010) return .p101_max;
    if (p010 >= p111 and p010 >= p000 and p010 >= p101) return .p010_max;
    return .equiprobable;
}

pub fn mostProbableSymbolOdd(p01: f64, p10: f64, bitwidth: u32) u64 {
    var mps: u64 = 0;
    const half = (bitwidth - 1) >> 1;
    switch (mostProbableTransitionPair(p01, p10)) {
        .p000_max, .equiprobable => { mps = 0; },
        .p111_max => {
            for (0..half) |_| { mps = (mps << 2) + 3; }
            mps = (mps << 1) + 1;
        },
        .p010_max => {
            for (0..half) |_| { mps = (mps << 2) + 1; }
            mps = (mps << 1) + 0;
        },
        .p101_max => {
            for (0..half) |_| { mps = (mps << 2) + 2; }
            mps = (mps << 1) + 1;
        },
    }
    return mps;
}

pub fn mostProbableSymbolEven(p01: f64, p10: f64, bitwidth: u32) u64 {
    var mps: u64 = 0;
    const p00 = 1.0 - p01;
    const p11 = 1.0 - p10;
    switch (mostProbableTransitionPair(p01, p10)) {
        .p000_max, .equiprobable => { mps = 0; },
        .p111_max => {
            for (0..@as(usize, bitwidth >> 1)) |_| { mps = (mps << 2) + 3; }
        },
        .p010_max => {
            if (bitwidth >= 2) {
                for (0..@as(usize, (bitwidth - 2) >> 1)) |_| { mps = (mps << 2) + 1; }
            }
            mps <<= 2;
            mps += if (p01 > p00) @as(u64, 1) else @as(u64, 0);
        },
        .p101_max => {
            if (bitwidth >= 2) {
                for (0..@as(usize, (bitwidth - 2) >> 1)) |_| { mps = (mps << 2) + 2; }
            }
            mps <<= 2;
            mps += if (p11 > p10) @as(u64, 3) else @as(u64, 2);
        },
    }
    return mps;
}

pub fn mostProbableSymbol(p01: f64, p10: f64, bitwidth: u32) u64 {
    return if ((bitwidth & 0x01) == 0x01)
        mostProbableSymbolOdd(p01, p10, bitwidth)
    else
        mostProbableSymbolEven(p01, p10, bitwidth);
}

pub const MaxProbResult = struct { probability: f64, mcv: u64 };

/// Compute the maximum symbol probability.
pub fn symbolMaxProbability(p01: f64, p10: f64, bitwidth: u32) MaxProbResult {
    const mu = p01 / (p10 + p01);
    const p0 = 1.0 - mu;
    const p1 = mu;
    const p00 = 1.0 - p01;
    const p11 = 1.0 - p10;
    const mps = mostProbableSymbol(p01, p10, bitwidth);

    var bits: [65]i32 = [_]i32{0} ** 65;
    bits[0] = 0;
    for (0..@as(usize, bitwidth)) |i| {
        const shift: u6 = @intCast(bitwidth - 1 - @as(u32, @intCast(i)));
        bits[i + 1] = @intCast((mps >> shift) & 0x01);
    }

    var p_0mps: f64 = 1.0;
    for (0..@as(usize, bitwidth)) |i| {
        if (bits[i] == 0 and bits[i + 1] == 0) { p_0mps *= p00; } else if (bits[i] == 0 and bits[i + 1] == 1) { p_0mps *= p01; } else if (bits[i] == 1 and bits[i + 1] == 0) { p_0mps *= p10; } else if (bits[i] == 1 and bits[i + 1] == 1) { p_0mps *= p11; }
    }

    bits[0] = 1;
    var p_1mps: f64 = 1.0;
    for (0..@as(usize, bitwidth)) |i| {
        if (bits[i] == 0 and bits[i + 1] == 0) { p_1mps *= p00; } else if (bits[i] == 0 and bits[i + 1] == 1) { p_1mps *= p01; } else if (bits[i] == 1 and bits[i + 1] == 0) { p_1mps *= p10; } else if (bits[i] == 1 and bits[i + 1] == 1) { p_1mps *= p11; }
    }

    return .{ .probability = (p0 * p_0mps) + (p1 * p_1mps), .mcv = mps };
}

pub const EntropyResult = struct { entropy: f64, mcv_prob: f64, mcv: u64 };

/// Convert Markov parameters to entropy per bit.
pub fn pToEntropy(p01: f64, p10: f64, bitwidth: u32) EntropyResult {
    const result = symbolMaxProbability(p01, p10, bitwidth);
    const ent = -math.log2(result.probability);
    return .{
        .entropy = ent / @as(f64, @floatFromInt(bitwidth)),
        .mcv_prob = result.probability,
        .mcv = result.mcv,
    };
}
