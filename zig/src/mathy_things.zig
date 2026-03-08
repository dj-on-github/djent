//! djent - A reimplementation of Fourmilab's ent with several improvements.
//!
//! Copyright (C) 2017  David Johnston
//! Translated to Zig 0.15.x, 2025
//!
//! Contact. David Johnston dj@deadhat.com

const std = @import("std");

/// Integer power function (wrapping on overflow, matching C behaviour).
pub fn ipow(base: u64, exp_val: u64) u64 {
    var result: u64 = 1;
    var b = base;
    var e = exp_val;
    while (e != 0) {
        if (e & 1 != 0) {
            result *%= b;
        }
        e >>= 1;
        b *%= b;
    }
    return result;
}

/// Standard normal CDF approximation.
pub fn zcdf(z: f64) f64 {
    if (z == 0.0) return 0.5;

    const y = @abs(z) / 2.0;

    if (y >= 3.0) return 0.0;

    var x: f64 = undefined;
    if (y < 1.0) {
        const w = y * y;
        x = 0.000124818987;
        x = x * w - 0.001075204047;
        x = x * w + 0.005198775019;
        x = x * w - 0.019198292004;
        x = x * w + 0.059054035642;
        x = x * w - 0.151968751364;
        x = x * w + 0.319152932694;
        x = x * w - 0.531923007300;
        x = x * w + 0.797884560593;
        x = x * 2.0 * y;
    } else {
        const yy = y - 2.0;
        x = -0.000045255659;
        x = x * yy + 0.000152529290;
        x = x * yy - 0.000019538132;
        x = x * yy - 0.000676904986;
        x = x * yy + 0.001390604284;
        x = x * yy - 0.000794620820;
        x = x * yy - 0.002034254874;
        x = x * yy + 0.006549791214;
        x = x * yy - 0.010557625006;
        x = x * yy + 0.011630447319;
        x = x * yy - 0.009279453341;
        x = x * yy + 0.005353579108;
        x = x * yy - 0.002141268741;
        x = x * yy + 0.000535310849;
        x = x * yy + 0.999936657524;
    }

    if (z > 0.0) {
        return (x / 2.0) + 0.5;
    } else {
        return 0.5 - (x / 2.0);
    }
}

const log_sqrt_pi: f64 = 0.5723649429247000870717135;
const i_sqrt_pi: f64 = 0.5641895835477562869480795;
const big_x: f64 = 20.0;

fn safeExp(x: f64) f64 {
    return if (x < -big_x) 0.0 else @exp(x);
}

/// Chi-square p-value.
pub fn chisqp(ax: f64, df: usize) f64 {
    const df_even = (df % 2) == 0;
    const x = ax;

    if (x <= 0.0 or df < 1) return 1.0;

    const a = x / 2.0;

    const y: f64 = if (df > 1) safeExp(-a) else 0.0;

    var s: f64 = if (df_even) y else 2.0 * zcdf(-@sqrt(x));

    if (df > 2) {
        const x_limit: f64 = @as(f64, @floatFromInt(df - 1)) / 2.0;
        var z: f64 = if (df_even) 1.0 else 0.5;

        if (a > big_x) {
            var e: f64 = if (df_even) 0.0 else log_sqrt_pi;
            const c = @log(a);
            while (z <= x_limit) {
                e = @log(z) + e;
                s += safeExp(c * z - a - e);
                z += 1.0;
            }
            return s;
        } else {
            var e: f64 = if (df_even) 1.0 else (i_sqrt_pi / @sqrt(a));
            var c: f64 = 0.0;
            while (z <= x_limit) {
                e = e * (a / z);
                c = c + e;
                z += 1.0;
            }
            return c * y + s;
        }
    } else {
        return s;
    }
}
