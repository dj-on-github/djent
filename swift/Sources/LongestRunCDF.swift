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

/// Return the probability of the longest run of heads being <= n
/// in a sequence of r uniform coin tosses.
///
/// The original C version used MPFR for arbitrary precision.
/// This Swift version uses Double, which is sufficient for typical
/// entropy analysis data sizes. For very large r values, a
/// multi-precision library (e.g. wrapping GMP/MPFR via C interop)
/// would be needed for full accuracy.
///
/// Formula:
///   first  = exp( -(r+1) / (2^(n+1) - n - 2) )
///   second = (2^(n+1) - 1) / (2^(n+1) - (n+2)/2)
///   result = first * second
func longestRunCDF(n: UInt, r: UInt) -> Double {
    let nd = Double(n)
    let rd = Double(r)

    // 2^(n+1)
    let twoToNPlus1 = pow(2.0, nd + 1.0)

    // first = exp( -(r+1) / (2^(n+1) - n - 2) )
    let topa = rd + 1.0
    let bottoma = twoToNPlus1 - nd - 2.0

    // Guard against division by zero for very small n
    guard bottoma != 0.0 else { return 0.0 }

    let first = exp(-topa / bottoma)

    // second = (2^(n+1) - 1) / (2^(n+1) - (n+2)/2)
    let topb = twoToNPlus1 - 1.0
    let nPlus2Over2 = (nd + 2.0) / 2.0
    let bottomb = twoToNPlus1 - nPlus2Over2

    guard bottomb != 0.0 else { return 0.0 }

    let second = topb / bottomb

    return first * second
}
