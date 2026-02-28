/*
    djent / djrandom - Markov 2-parameter model.
    
    Copyright (C) 2017  David Johnston
    Translated to Swift 2026

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Contact. David Johnston dj@deadhat.com
*/

import Foundation

// A library for converting between points, SCC, bias and entropy
// with the 2-parameter Markov model.

enum TransitionPairType {
    case equiprobable
    case p000Max
    case p111Max
    case p101Max
    case p010Max
}

var markovVerboseMode = 0

/// Format a symbol as a binary string of given bit width
func symbolText(_ x: UInt64, bitwidth: Int) -> String {
    var s = ""
    for i in 0..<bitwidth {
        let bit = (x >> (bitwidth - 1 - i)) & 0x01
        s += bit == 0 ? "0" : "1"
    }
    return s
}

/// Make two probability density functions for all 2^bitwidth symbols.
/// One for when the previous bit is 0, one for when it is 1.
func makePDF(p01: Double, p10: Double, bitwidth: Int) -> (table0: [Double], table1: [Double]) {
    let size = 1 << bitwidth
    let p00 = 1.0 - p01
    let p11 = 1.0 - p10
    var table0 = [Double](repeating: 0.0, count: size)
    var table1 = [Double](repeating: 0.0, count: size)
    var sum0 = 0.0
    var sum1 = 0.0

    for x in 0..<size {
        if p01 == 0.5 && p10 == 0.5 {
            let uniform = 1.0 / Double(size)
            table0[x] = uniform
            table1[x] = uniform
            sum0 += uniform
            sum1 = sum0
        } else {
            var plist0 = 1.0
            var plist1 = 1.0

            // First bit with previous last bit
            if (x & 0x1) == 0 {
                plist0 *= p00
                plist1 *= p10
            } else {
                plist0 *= p01
                plist1 *= p11
            }

            for i in 0..<(bitwidth - 1) {
                let bp = (x >> i) & 0x3
                switch bp {
                case 0:
                    plist0 *= p00; plist1 *= p00
                case 1:
                    plist0 *= p10; plist1 *= p10
                case 2:
                    plist0 *= p01; plist1 *= p01
                case 3:
                    plist0 *= p11; plist1 *= p11
                default:
                    break
                }
            }

            table0[x] = plist0
            table1[x] = plist1
            sum0 += plist0
            sum1 += plist1
        }
    }

    // Normalize
    for i in 0..<256 {
        guard i < size else { break }
        table0[i] /= sum0
        table1[i] /= sum1
    }

    return (table0, table1)
}

/// Make two cumulative density functions for all 2^bitwidth symbols.
func makeCDF(p01: Double, p10: Double, bitwidth: Int) -> (table0: [Double], table1: [Double]) {
    let size = 1 << bitwidth
    let p00 = 1.0 - p01
    let p11 = 1.0 - p10
    var table0 = [Double](repeating: 0.0, count: size)
    var table1 = [Double](repeating: 0.0, count: size)

    for x in 0..<size {
        if p01 == 0.5 && p10 == 0.5 {
            let uniform = 1.0 / Double(size)
            if x == 0 {
                table0[x] = uniform
                table1[x] = uniform
            } else {
                table0[x] = table0[x - 1] + uniform
                table1[x] = table1[x - 1] + uniform
            }
        } else {
            var plist0 = 1.0
            var plist1 = 1.0

            if (x & 0x1) == 0 {
                plist0 *= p00
                plist1 *= p10
            } else {
                plist0 *= p01
                plist1 *= p11
            }

            for i in 0..<(bitwidth - 1) {
                let bp = (x >> i) & 0x3
                switch bp {
                case 0:
                    plist0 *= p00; plist1 *= p00
                case 1:
                    plist0 *= p10; plist1 *= p10
                case 2:
                    plist0 *= p01; plist1 *= p01
                case 3:
                    plist0 *= p11; plist1 *= p11
                default:
                    break
                }
            }

            if x == 0 {
                table0[x] = plist0
                table1[x] = plist1
            } else {
                table0[x] = table0[x - 1] + plist0
                table1[x] = table1[x - 1] + plist1
            }
        }
    }

    // Normalize
    let max0 = table0[size - 1]
    let max1 = table1[size - 1]
    for i in 0..<size {
        table0[i] /= max0
        table1[i] /= max1
    }

    return (table0, table1)
}

/// Build sampling tables for Markov generation
func makeSampleTable(p01: Double, p10: Double, bitwidth: Int) -> (sampleTable0: [Int], sampleTable1: [Int]) {
    let (cdfTable0, cdfTable1) = makeCDF(p01: p01, p10: p10, bitwidth: bitwidth)
    let (_, _) = makePDF(p01: p01, p10: p10, bitwidth: bitwidth)

    let tableSize = 1 << 20
    var st0 = [Int](repeating: 0, count: tableSize)
    var st1 = [Int](repeating: 0, count: tableSize)

    // Populate the 1M table with symbols according to the CDF
    var baseIndex = 0
    for x in 0..<256 {
        let floatPos = cdfTable0[x]
        let index = Int(floatPos * Double(tableSize))
        for i in baseIndex..<index {
            st0[i] = x
        }
        baseIndex = index
    }

    baseIndex = 0
    for x in 0..<256 {
        let floatPos = cdfTable1[x]
        let index = Int(floatPos * Double(tableSize))
        for i in baseIndex..<index {
            st1[i] = x
        }
        baseIndex = index
    }

    return (st0, st1)
}

/// Compute the symbol probability for the Markov 2-parameter model
func symbolProb(p01: Double, p10: Double, x: UInt64, bitwidth: Int) -> Double {
    let p00 = 1.0 - p01
    let p11 = 1.0 - p10
    let mu = p01 / (p10 + p01)
    let p0 = 1.0 - mu
    let p1 = mu

    if p01 == 0.5 && p10 == 0.5 { return 1.0 }

    var plist0 = 1.0
    var plist1 = 1.0

    if ((x >> (bitwidth - 1)) & 0x1) == 0 {
        plist0 *= p00
        plist1 *= p10
    } else {
        plist0 *= p01
        plist1 *= p11
    }

    for i in 0..<(bitwidth - 2) {
        let bp = Int((x >> (bitwidth - 2 - i)) & 0x3)
        switch bp {
        case 0: plist0 *= p00; plist1 *= p00
        case 1: plist0 *= p01; plist1 *= p01
        case 2: plist0 *= p10; plist1 *= p10
        case 3: plist0 *= p11; plist1 *= p11
        default: break
        }
    }

    return (p0 * plist0) + (p1 * plist1)
}

func maxDouble(_ x: Double, _ y: Double) -> Double {
    return Swift.max(x, y)
}

func mkSymbol(prefix: Int, tbp: Int, postfix: Int, bitwidth: Int) -> UInt64 {
    let rep = (bitwidth - 2) / 2
    var pattern = UInt64(prefix)
    for _ in 0..<rep {
        pattern = (pattern << 2) + UInt64(tbp)
    }
    pattern = (pattern << 1) + UInt64(postfix)
    return pattern
}

func mkSymbolNoPostfix(prefix: Int, tbp: Int, bitwidth: Int) -> UInt64 {
    var pattern = UInt64(prefix)
    for _ in 0..<((bitwidth - 1) / 2) {
        pattern = (pattern << 2) + UInt64(tbp)
    }
    return pattern
}

func mostProbableTransitionPair(p01: Double, p10: Double) -> TransitionPairType {
    let mu = p01 / (p10 + p01)
    let p0 = 1.0 - mu
    let p1 = mu
    let p00 = 1.0 - p01
    let p11 = 1.0 - p10

    let p010 = p0 * p01 * p10
    let p101 = p1 * p10 * p01
    let p000 = p0 * p00 * p00
    let p111 = p1 * p11 * p11

    if p111 >= p000 && p111 >= p101 && p111 >= p010 {
        return .p111Max
    } else if p000 >= p111 && p000 >= p101 && p000 >= p010 {
        return .p000Max
    } else if p101 >= p111 && p101 >= p000 && p101 >= p010 {
        return .p101Max
    } else if p010 >= p111 && p010 >= p000 && p010 >= p101 {
        return .p010Max
    }
    return .equiprobable
}

func mostProbableSymbolOdd(p01: Double, p10: Double, bitwidth: Int) -> UInt64 {
    var mps: UInt64 = 0
    let tp = mostProbableTransitionPair(p01: p01, p10: p10)

    switch tp {
    case .p000Max:
        mps = 0
    case .p111Max:
        for _ in 0..<((bitwidth - 1) >> 1) {
            mps = (mps << 2) + 3
        }
        mps = (mps << 1) + 1
    case .p010Max:
        for _ in 0..<((bitwidth - 1) >> 1) {
            mps = (mps << 2) + 1
        }
        mps = (mps << 1) + 0
    case .p101Max:
        for _ in 0..<((bitwidth - 1) >> 1) {
            mps = (mps << 2) + 2
        }
        mps = (mps << 1) + 1
    case .equiprobable:
        mps = 0
    }
    return mps
}

func mostProbableSymbolEven(p01: Double, p10: Double, bitwidth: Int) -> UInt64 {
    var mps: UInt64 = 0
    let p00 = 1.0 - p01
    let p11 = 1.0 - p10
    let tp = mostProbableTransitionPair(p01: p01, p10: p10)

    switch tp {
    case .p000Max:
        mps = 0
    case .p111Max:
        for _ in 0..<(bitwidth >> 1) {
            mps = (mps << 2) + 3
        }
    case .p010Max:
        for _ in 0..<((bitwidth - 2) >> 1) {
            mps = (mps << 2) + 1
        }
        mps <<= 2
        mps += p01 > p00 ? 1 : 0
    case .p101Max:
        for _ in 0..<((bitwidth - 2) >> 1) {
            mps = (mps << 2) + 2
        }
        mps <<= 2
        mps += p11 > p10 ? 3 : 2
    case .equiprobable:
        mps = 0
    }
    return mps
}

func mostProbableSymbol(p01: Double, p10: Double, bitwidth: Int) -> UInt64 {
    let mps: UInt64
    if (bitwidth & 0x01) == 0x01 {
        mps = mostProbableSymbolOdd(p01: p01, p10: p10, bitwidth: bitwidth)
    } else {
        mps = mostProbableSymbolEven(p01: p01, p10: p10, bitwidth: bitwidth)
    }

    if markovVerboseMode > 1 {
        fputs("   MCV = 0x\(String(mps, radix: 16))\n", stderr)
    }
    return mps
}

/// Compute the maximum symbol probability for the given Markov model parameters
func symbolMaxProbability(p01: Double, p10: Double, bitwidth: Int) -> (probability: Double, mcv: UInt64) {
    let mu = p01 / (p10 + p01)
    let p0 = 1.0 - mu
    let p1 = mu
    let p00 = 1.0 - p01
    let p11 = 1.0 - p10

    let mps = mostProbableSymbol(p01: p01, p10: p10, bitwidth: bitwidth)

    // Unpack symbol bits into an array
    var bits = [Int](repeating: 0, count: bitwidth + 1)

    // First with x[-1] = 0
    bits[0] = 0
    for i in 0..<bitwidth {
        bits[i + 1] = Int((mps >> (bitwidth - 1 - i)) & 0x01)
    }

    if markovVerboseMode > 1 {
        let bitStr = bits.prefix(bitwidth + 1).map { String($0) }.joined()
        fputs("   unrolled bits 0 prefix = \(bitStr)\n", stderr)
    }

    // Compute probability with prefix 0
    var p0mps = 1.0
    for i in 0..<bitwidth {
        if bits[i] == 0 && bits[i + 1] == 0      { p0mps *= p00 }
        else if bits[i] == 0 && bits[i + 1] == 1  { p0mps *= p01 }
        else if bits[i] == 1 && bits[i + 1] == 0  { p0mps *= p10 }
        else if bits[i] == 1 && bits[i + 1] == 1  { p0mps *= p11 }
    }

    // Then with x[-1] = 1
    bits[0] = 1

    if markovVerboseMode > 1 {
        let bitStr = bits.prefix(bitwidth + 1).map { String($0) }.joined()
        fputs("   unrolled bits 1 prefix = \(bitStr)\n", stderr)
    }

    var p1mps = 1.0
    for i in 0..<bitwidth {
        if bits[i] == 0 && bits[i + 1] == 0      { p1mps *= p00 }
        else if bits[i] == 0 && bits[i + 1] == 1  { p1mps *= p01 }
        else if bits[i] == 1 && bits[i + 1] == 0  { p1mps *= p10 }
        else if bits[i] == 1 && bits[i + 1] == 1  { p1mps *= p11 }
    }

    let pMps = (p0 * p0mps) + (p1 * p1mps)
    return (pMps, mps)
}

/// Convert Markov parameters to entropy per bit
func pToEntropy(p01: Double, p10: Double, bitwidth: Int) -> (entropy: Double, mcvProb: Double, mcv: UInt64) {
    let (smp, mcv) = symbolMaxProbability(p01: p01, p10: p10, bitwidth: bitwidth)
    let ent = -log2(smp)
    return (ent / Double(bitwidth), smp, mcv)
}

func near(_ x: Double, _ y: Double, epsilon: Double) -> Bool {
    return y > x - epsilon && y < x + epsilon
}
