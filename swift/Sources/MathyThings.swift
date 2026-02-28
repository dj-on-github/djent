/*
    djent - A reimplementation of Fourmilab's ent with several improvements.
    
    Copyright (C) 2017  David Johnston
    Translated to Swift 2026

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Contact. David Johnston dj@deadhat.com
*/

import Foundation

/// Integer power function
func ipow(_ base: UInt64, _ exp: UInt64) -> UInt64 {
    var result: UInt64 = 1
    var base = base
    var exp = exp
    while exp != 0 {
        if exp & 1 != 0 {
            result &*= base
        }
        exp >>= 1
        base &*= base
    }
    return result
}

// MARK: - Chi-Square P-value computation

/// Standard normal cumulative distribution function approximation
func zcdf(_ z: Double) -> Double {
    if z == 0.0 { return 0.5 }

    let y = abs(z) / 2.0

    if y >= 3.0 { return 0.0 }

    let x: Double
    if y < 1.0 {
        let w = y * y
        var v = 0.000124818987
        v = v * w - 0.001075204047
        v = v * w + 0.005198775019
        v = v * w - 0.019198292004
        v = v * w + 0.059054035642
        v = v * w - 0.151968751364
        v = v * w + 0.319152932694
        v = v * w - 0.531923007300
        v = v * w + 0.797884560593
        x = v * 2.0 * y
    } else {
        let yy = y - 2.0
        var v = -0.000045255659
        v = v * yy + 0.000152529290
        v = v * yy - 0.000019538132
        v = v * yy - 0.000676904986
        v = v * yy + 0.001390604284
        v = v * yy - 0.000794620820
        v = v * yy - 0.002034254874
        v = v * yy + 0.006549791214
        v = v * yy - 0.010557625006
        v = v * yy + 0.011630447319
        v = v * yy - 0.009279453341
        v = v * yy + 0.005353579108
        v = v * yy - 0.002141268741
        v = v * yy + 0.000535310849
        v = v * yy + 0.999936657524
        x = v
    }

    if z > 0.0 {
        return (x / 2.0) + 0.5
    } else {
        return 0.5 - (x / 2.0)
    }
}

private let logSqrtPi = 0.5723649429247000870717135  // log(sqrt(pi))
private let iSqrtPi   = 0.5641895835477562869480795  // 1/sqrt(pi)
private let bigX      = 20.0                          // max value for exp(x)

private func safeExp(_ x: Double) -> Double {
    return x < -bigX ? 0.0 : exp(x)
}

/// Chi-square p-value
func chisqp(_ ax: Double, _ df: Int) -> Double {
    let dfEven = (df % 2) == 0

    let x = ax
    if x <= 0.0 || df < 1 { return 1.0 }

    let a = x / 2.0

    let y: Double = df > 1 ? safeExp(-a) : 0.0

    var s: Double
    if dfEven {
        s = y
    } else {
        s = 2.0 * zcdf(-sqrt(x))
    }

    if df > 2 {
        let xLimit = Double(df - 1) / 2.0
        var z: Double = dfEven ? 1.0 : 0.5

        if a > bigX {
            var e: Double = dfEven ? 0.0 : logSqrtPi
            let c = log(a)
            while z <= xLimit {
                e = log(z) + e
                s += safeExp(c * z - a - e)
                z += 1.0
            }
            return s
        } else {
            var e: Double = dfEven ? 1.0 : (iSqrtPi / sqrt(a))
            var c = 0.0
            while z <= xLimit {
                e = e * (a / z)
                c = c + e
                z += 1.0
            }
            return c * y + s
        }
    } else {
        return s
    }
}
