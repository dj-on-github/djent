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

// MARK: - Constants

private let queueSize = 4096
private let buffSize  = 2048

// MARK: - Byte Reverse Table

private let byteReverseTable: [UInt8] = [
    0x00,0x80,0x40,0xC0,0x20,0xA0,0x60,0xE0,0x10,0x90,0x50,0xD0,0x30,0xB0,0x70,0xF0,
    0x08,0x88,0x48,0xC8,0x28,0xA8,0x68,0xE8,0x18,0x98,0x58,0xD8,0x38,0xB8,0x78,0xF8,
    0x04,0x84,0x44,0xC4,0x24,0xA4,0x64,0xE4,0x14,0x94,0x54,0xD4,0x34,0xB4,0x74,0xF4,
    0x0C,0x8C,0x4C,0xCC,0x2C,0xAC,0x6C,0xEC,0x1C,0x9C,0x5C,0xDC,0x3C,0xBC,0x7C,0xFC,
    0x02,0x82,0x42,0xC2,0x22,0xA2,0x62,0xE2,0x12,0x92,0x52,0xD2,0x32,0xB2,0x72,0xF2,
    0x0A,0x8A,0x4A,0xCA,0x2A,0xAA,0x6A,0xEA,0x1A,0x9A,0x5A,0xDA,0x3A,0xBA,0x7A,0xFA,
    0x06,0x86,0x46,0xC6,0x26,0xA6,0x66,0xE6,0x16,0x96,0x56,0xD6,0x36,0xB6,0x76,0xF6,
    0x0E,0x8E,0x4E,0xCE,0x2E,0xAE,0x6E,0xEE,0x1E,0x9E,0x5E,0xDE,0x3E,0xBE,0x7E,0xFE,
    0x01,0x81,0x41,0xC1,0x21,0xA1,0x61,0xE1,0x11,0x91,0x51,0xD1,0x31,0xB1,0x71,0xF1,
    0x09,0x89,0x49,0xC9,0x29,0xA9,0x69,0xE9,0x19,0x99,0x59,0xD9,0x39,0xB9,0x79,0xF9,
    0x05,0x85,0x45,0xC5,0x25,0xA5,0x65,0xE5,0x15,0x95,0x55,0xD5,0x35,0xB5,0x75,0xF5,
    0x0D,0x8D,0x4D,0xCD,0x2D,0xAD,0x6D,0xED,0x1D,0x9D,0x5D,0xDD,0x3D,0xBD,0x7D,0xFD,
    0x03,0x83,0x43,0xC3,0x23,0xA3,0x63,0xE3,0x13,0x93,0x53,0xD3,0x33,0xB3,0x73,0xF3,
    0x0B,0x8B,0x4B,0xCB,0x2B,0xAB,0x6B,0xEB,0x1B,0x9B,0x5B,0xDB,0x3B,0xBB,0x7B,0xFB,
    0x07,0x87,0x47,0xC7,0x27,0xA7,0x67,0xE7,0x17,0x97,0x57,0xD7,0x37,0xB7,0x77,0xF7,
    0x0F,0x8F,0x4F,0xCF,0x2F,0xAF,0x6F,0xEF,0x1F,0x9F,0x5F,0xDF,0x3F,0xBF,0x7F,0xFF,
]

// MARK: - Command-Line Options

struct Options {
    var symbolLength: UInt = 8
    var hexMode: Bool = true
    var printOccurrence: Bool = false
    var printLongest: Bool = false
    var fold: Bool = false
    var terse: Bool = false
    var useStdin: Bool = true
    var sccWrap: Bool = false
    var lagN: Int = 1
    var usingInputListFile: Bool = false
    var inputListFilename: String = ""
    var suppressHeader: Bool = false
    var byteReverse: Bool = false
    var parseFilename: Bool = false
    var wordReverse: Bool = false
    var entExact: Bool = false
    var gotSkip: Bool = false
    var skipAmount: Int = 0
    var gotSubstring: Bool = false
    var substring: Int = 0
}

// MARK: - Analysis Results

struct AnalysisResults {
    var mean: Double = 0.0
    var chisqCount: UInt64 = 0
    var chisqDistribution: Double = 0.0
    var chisqPercent: Double = 0.0
    var entropy: Double = 0.0
    var minEntropy: Double = 0.0
    var minEntropySymbol: UInt32 = 0
    var pi: Double = 0.0
    var piErr: Double = 0.0
    var compression: Double = 0.0
    var scc: Double = 0.0
    var p01: Double = 0.0
    var p10: Double = 0.0
    var longestPValue: Double = 0.0
}

// MARK: - Byte Queue (FIFO)

/// A FIFO queue that accepts bytes and dispenses symbols of arbitrary bit width.
class ByteQueue {
    private var queue = [UInt8](repeating: 0, count: queueSize)
    private var queueStart: Int = 0
    private var queueEnd: Int = 0
    var queueUsed: Int = 0

    private var currentByte: UInt = 0
    private var bitsUsedFromByte: UInt = 0
    private var gotByte: Bool = false
    private var currentSymbol: Int64 = 0
    private var bitsInCurrentSymbol: UInt = 0
    var symbolMask: UInt64 = 0

    func initialize(symbolLength: UInt) {
        queueStart = 0
        queueEnd = 0
        queueUsed = 0
        gotByte = false
        currentByte = 0
        bitsUsedFromByte = 0
        currentSymbol = 0
        bitsInCurrentSymbol = 0
        for i in 0..<queueSize { queue[i] = 0 }
        symbolMask = ipow(2, UInt64(symbolLength)) - 1
    }

    /// Push a byte into the queue, optionally bit-reversing it
    func push(_ byte: UInt8, byteReverse: Bool) {
        let b = byteReverse ? byteReverseTable[Int(byte)] : byte
        queue[queueEnd] = b
        queueEnd = (queueEnd + 1) % queueSize
        queueUsed += 1
    }

    /// Space remaining in the queue
    var spaceRemaining: Int {
        return queueSize - queueUsed
    }

    /// Pull a symbol of the given bit width from the queue. Returns -1 if insufficient data.
    func getSymbol(_ symbolLength: UInt) -> Int64 {
        currentSymbol = 0

        // Get a byte if we don't have one
        if !gotByte {
            if UInt(queueUsed) * 8 < symbolLength { return -1 }
            currentByte = UInt(queue[queueStart])
            queueStart = (queueStart + 1) % queueSize
            queueUsed -= 1
            bitsUsedFromByte = 0
            gotByte = true
        }

        // Optimized: single bit
        if symbolLength == 1 {
            currentSymbol = Int64((currentByte & 0x80) >> 7)
            currentByte = (currentByte << 1) & 0xFF
            bitsUsedFromByte += 1
            if bitsUsedFromByte == 8 {
                gotByte = false
                bitsUsedFromByte = 0
            }
            return currentSymbol
        }

        // Optimized: full byte
        if symbolLength == 8 {
            currentSymbol = Int64(currentByte)
            gotByte = false
            return currentSymbol
        }

        // General case: bit-by-bit extraction
        bitsInCurrentSymbol = 0
        repeat {
            let temp = (currentByte & 0x80) >> 7
            currentByte = (currentByte << 1) & 0xFF
            bitsUsedFromByte += 1
            if bitsUsedFromByte == 8 {
                gotByte = false
                bitsUsedFromByte = 0
            }
            currentSymbol = Int64((UInt64(currentSymbol << 1) | UInt64(temp)) & symbolMask)
            bitsInCurrentSymbol += 1

            // Fetch a new byte if needed
            if !gotByte {
                if UInt(queueUsed) * 8 < symbolLength { return -1 }
                currentByte = UInt(queue[queueStart])
                queueStart = (queueStart + 1) % queueSize
                queueUsed -= 1
                bitsUsedFromByte = 0
                gotByte = true
            }
        } while bitsInCurrentSymbol < symbolLength

        return currentSymbol
    }
}

// MARK: - Hex-to-Binary Converter

class HexConverter {
    private var hexpair: [UInt8] = [0, 0]
    private var state: Int = 0

    func reset() {
        hexpair = [0, 0]
        state = 0
    }

    private func isHex(_ c: UInt8) -> Bool {
        return (c >= 0x30 && c <= 0x39) ||  // 0-9
               (c >= 0x41 && c <= 0x46) ||  // A-F
               (c >= 0x61 && c <= 0x66)     // a-f
    }

    private func isHexOrX(_ c: UInt8) -> Bool {
        return isHex(c) || c == UInt8(ascii: "x")
    }

    /// Convert hex text in buffer to binary in-place. Returns the number of output bytes.
    func convert(_ buffer: inout [UInt8], length: Int) -> Int {
        var outPos = 0
        var scanPos = 0

        while scanPos < length {
            if state == 0 {
                let c = buffer[scanPos]
                if isHex(c) {
                    hexpair[0] = c
                    state = 1
                }
                scanPos += 1
            } else if state == 1 {
                let c = buffer[scanPos]
                if (isHexOrX(c) && hexpair[0] == UInt8(ascii: "0")) || isHex(c) {
                    hexpair[1] = c
                    state = 2
                }
                scanPos += 1
            } else if state == 2 {
                if hexpair[0] == UInt8(ascii: "0") && hexpair[1] == UInt8(ascii: "x") {
                    state = 0
                } else {
                    // Valid hex pair
                    var nybble = 0
                    let h0 = Int(hexpair[0])
                    if h0 > 47 && h0 < 58       { nybble = h0 - 48 }
                    else if h0 > 64 && h0 < 71   { nybble = h0 - 55 }
                    else if h0 > 96 && h0 < 103   { nybble = h0 - 87 }
                    nybble <<= 4

                    var byte = 0
                    let h1 = Int(hexpair[1])
                    if h1 > 47 && h1 < 58        { byte = nybble + h1 - 48 }
                    else if h1 > 64 && h1 < 71    { byte = nybble + h1 - 55 }
                    else if h1 > 96 && h1 < 103    { byte = nybble + h1 - 87 }

                    buffer[outPos] = UInt8(byte)
                    outPos += 1
                    state = 0
                }
            }
        }
        return outPos
    }
}

// MARK: - Entropy Analyzer

class EntropyAnalyzer {
    let opts: Options
    let queue = ByteQueue()
    let hexConverter = HexConverter()

    // Analysis state
    var symbolCount: UInt64 = 0
    var meanTotal: UInt64 = 0
    var filebytes: UInt64 = 0
    var ent: Double = 0.0

    // Occurrence tracking
    var occurrenceSize: UInt64 = 0
    var occurrenceCount: [UInt64] = []
    var occurrenceTotal: UInt64 = 0
    var noOccurrenceSpace: Bool = false

    // Longest run tracking
    var longestSize: UInt64 = 0
    var longestCount: [UInt64] = []
    var longestTotal: UInt64 = 0
    var noLongestSpace: Bool = false
    var longestLastSymbol: UInt64 = 0
    var longestRun: UInt64 = 0
    var longestLongest: UInt64 = 0
    var longestLongestSymbol: UInt64 = 0
    var longestPosition: UInt64 = 0
    var longestNewPos: UInt64 = 0
    var longestBytePos: UInt64 = 0

    // Chi-square
    var chisq: Double = 0.0
    var chisqSum: Double = 0.0
    var chisqProb: [Double] = []

    // Monte Carlo Pi
    var mp: UInt64 = 0
    var montyTotalCount: UInt64 = 0
    var montyInsideCount: UInt64 = 0
    var radiusSquared: Double = 0.0
    var positionX: Double = 0.0
    var positionY: Double = 0.0
    var monte: [UInt64] = [0, 0, 0, 0, 0, 0]

    // Serial Correlation
    var t1: UInt64 = 0
    var t2: UInt64 = 0
    var t3: UInt64 = 0
    var sccFirst: Bool = true
    var firstSymbol: UInt64 = 0
    var sccPrevious: UInt64 = 0
    var sccCount: UInt64 = 0
    var sccFIFO = [UInt64](repeating: 0, count: 256)
    var sccFirstLagN = [UInt64](repeating: 0, count: 256)
    var aeqbCount: UInt64 = 0
    var meanCount: UInt64 = 0
    var otherSCC: Double = 0.0

    // Results
    var results = AnalysisResults()
    var markovEntropy: Double = 0.0

    // Word-reverse buffer
    var buffer2 = [UInt8](repeating: 0, count: buffSize + 4)
    var buffer2Size: Int = 0

    init(options: Options) {
        self.opts = options
    }

    // MARK: - Initialization

    func initAll() {
        symbolCount = 0
        queue.initialize(symbolLength: opts.symbolLength)
        hexConverter.reset()
        initMean()
        initEntropy()
        initOccurrences()
        initLongest()
        initChisq()
        initMonteCarlo()
        initSCC()
        filebytes = 0
        buffer2Size = 0
    }

    func initMean() {
        meanTotal = 0
    }

    func initEntropy() {
        ent = 0.0
    }

    func initOccurrences() {
        noOccurrenceSpace = false
        occurrenceTotal = 0
        if opts.symbolLength > 32 {
            fputs("Error, symbol length cannot be longer than 32 bits for occurrence count table\n", stderr)
            exit(1)
        }
        occurrenceSize = ipow(2, UInt64(opts.symbolLength))
        occurrenceCount = [UInt64](repeating: 0, count: Int(occurrenceSize))
    }

    func initLongest() {
        noLongestSpace = false
        longestTotal = 0
        longestPosition = 0
        if opts.symbolLength > 32 {
            fputs("Error, symbol length cannot be longer than 32 bits for longest count table\n", stderr)
            exit(1)
        }
        longestSize = ipow(2, UInt64(opts.symbolLength))
        longestCount = [UInt64](repeating: 0, count: Int(longestSize))
        longestLastSymbol = 0
        longestRun = 0
        longestLongest = 0
        longestLongestSymbol = 0
        longestNewPos = 0
        longestBytePos = 0
    }

    func initChisq() {
        chisq = 0.0
        chisqSum = 0.0
        chisqProb = [Double](repeating: 0.0, count: Int(occurrenceSize))
    }

    func initMonteCarlo() {
        mp = 0
        montyTotalCount = 0
        montyInsideCount = 0
        radiusSquared = (256.0 * 256.0 * 256.0) - 1.0
        radiusSquared = radiusSquared * radiusSquared
        monte = [0, 0, 0, 0, 0, 0]
    }

    func initSCC() {
        t1 = 0; t2 = 0; t3 = 0
        sccFirst = true
        sccPrevious = 0
        sccCount = 0
        firstSymbol = 0
        sccFIFO = [UInt64](repeating: 0, count: 256)
        sccFirstLagN = [UInt64](repeating: 0, count: 256)
        aeqbCount = 0
        meanCount = 0
    }

    // MARK: - Update routines

    func updateMean(_ symbol: UInt64) {
        meanTotal += symbol
    }

    func updateOccurrences(_ symbol: UInt64) {
        occurrenceCount[Int(symbol)] += 1
        occurrenceTotal += 1
    }

    func updateLongest(_ symbol: UInt64, symbolPos: UInt64) {
        let symbolBytePos = (UInt64(opts.symbolLength) * symbolPos) / 8

        if symbol == longestLastSymbol {
            longestRun += 1
            if longestRun > longestCount[Int(symbol)] {
                longestCount[Int(symbol)] = longestRun
            }
            if longestRun > longestLongest {
                longestLongest = longestRun
                longestLongestSymbol = symbol
                longestPosition = longestNewPos
            }
        } else {
            longestRun = 1
            longestLastSymbol = symbol
            longestNewPos = symbolBytePos
        }
    }

    func updateMonteCarlo(_ symbol: UInt8) {
        monte[Int(mp)] = UInt64(symbol)
        mp += 1

        if mp > 5 {
            mp = 0
            montyTotalCount += 1
            positionX = 0.0
            positionY = 0.0
            for mj in 0..<3 {
                positionX = (positionX * 256.0) + Double(monte[mj])
                positionY = (positionY * 256.0) + Double(monte[3 + mj])
            }
            if (positionX * positionX) + (positionY * positionY) <= radiusSquared {
                montyInsideCount += 1
            }
        }
    }

    func updateSCC(_ symbol: UInt64) {
        let lagn = opts.lagN
        if lagn == 1 {
            sccCount += 1
            if sccFirst {
                sccFirst = false
                firstSymbol = symbol
            } else {
                t1 += sccPrevious * symbol
                if sccPrevious == symbol { aeqbCount += 1 }
            }
            meanCount += symbol
            t2 += symbol * symbol
            t3 += symbol
            sccPrevious = symbol
        } else {
            sccCount += 1
            if sccCount <= UInt64(lagn) {
                sccFIFO[Int(sccCount - 1)] = symbol
                sccFirstLagN[Int(sccCount - 1)] = symbol
            } else {
                t1 += sccFIFO[0] * symbol
                if sccFIFO[0] == symbol { aeqbCount += 1 }
                for i in 0..<lagn {
                    sccFIFO[i] = sccFIFO[i + 1]
                }
                meanCount += symbol
                sccFIFO[lagn] = symbol
                t2 += symbol * symbol
                t3 += symbol
            }
        }
    }

    // MARK: - Finalization routines

    func finalizeMean() {
        results.mean = Double(meanTotal) / Double(symbolCount)
    }

    func finalizeEntropy() {
        ent = 0.0
        for i in 0..<Int(occurrenceSize) {
            if chisqProb[i] > 0.0 {
                ent += chisqProb[i] * log10(1.0 / chisqProb[i]) * 3.32192809488736234787
            }
        }
        results.entropy = ent
    }

    func finalizeOccurrences() {
        var maxC: UInt64 = 0
        var maxSymbol: UInt32 = 0
        for i in 0..<Int(occurrenceSize) {
            if occurrenceCount[i] > maxC {
                maxC = occurrenceCount[i]
                maxSymbol = UInt32(i)
            }
        }
        let maxP = Double(maxC) / Double(occurrenceTotal)
        let maxPEnt = (-log10(maxP) / log10(2.0)) / Double(opts.symbolLength)
        results.minEntropy = maxPEnt
        results.minEntropySymbol = maxSymbol

        if !opts.terse {
            print("   Min Entropy (by max occurrence of symbol \(String(maxSymbol, radix: 16))) = \(maxPEnt)")
        }
    }

    func finalizeLongest() {
        results.longestPValue = longestRunCDF(n: UInt(longestLongest), r: UInt(symbolCount))

        if opts.symbolLength != 8 {
            longestBytePos = (symbolCount * UInt64(opts.symbolLength)) / 8
        } else {
            longestBytePos = symbolCount
        }
    }

    func finalizeChisq() {
        let expected = Double(occurrenceTotal) / Double(occurrenceSize)
        chisq = 0.0
        chisqSum = 0.0
        for i in 0..<Int(occurrenceSize) {
            let diff = Double(occurrenceCount[i]) - expected
            chisqProb[i] = Double(occurrenceCount[i]) / Double(occurrenceTotal)
            chisq += (diff * diff) / expected
            chisqSum += Double(UInt64(i) * occurrenceCount[i])
        }
        let finalProb = chisqp(chisq, Int(occurrenceSize) - 1)
        results.chisqCount = occurrenceTotal
        results.chisqDistribution = chisq
        results.chisqPercent = finalProb * 100.0
    }

    func finalizeMonteCarlo() {
        let montepi = 4.0 * (Double(montyInsideCount) / Double(montyTotalCount))
        let pierr = (abs(Double.pi - montepi) / Double.pi) * 100.0
        results.pi = montepi
        results.piErr = pierr
    }

    func finalizeCompression() {
        let sl = Double(opts.symbolLength)
        results.compression = (100.0 * (sl - ent)) / sl
    }

    func finalizeSCC() {
        let lagn = opts.lagN

        if opts.sccWrap {
            if lagn == 1 {
                t1 += sccPrevious * firstSymbol
                t2 += firstSymbol * firstSymbol
                t3 += firstSymbol
            } else {
                for i in 0..<lagn {
                    t1 += sccFIFO[i] * sccFirstLagN[i]
                    t2 += sccFirstLagN[i] * sccFirstLagN[i]
                    t3 += sccFirstLagN[i]
                }
            }
        } else {
            sccCount -= UInt64(lagn)
        }

        let top = Int64(sccCount &* t1) &- Int64(t3 &* t3)
        let bottom = Int64(sccCount &* t2) &- Int64(t3 &* t3)
        let scc = Double(top) / Double(bottom)

        results.scc = scc

        // A=B count based SCC estimate
        let bias = Double(t3) / Double(sccCount)
        let paeqb = Double(aeqbCount) / Double(sccCount)
        otherSCC = (2.0 * paeqb) - 1.0
        otherSCC = otherSCC * pow(1.0 - (2.0 * abs(bias - 0.5)), 2)
    }

    func computeMarkov() {
        let p01 = results.mean * (1.0 - results.scc)
        let p10 = (1.0 - results.mean) * (1.0 - results.scc)
        results.p01 = p01
        results.p10 = p10

        let (entropy, _, _) = pToEntropy(p01: p01, p10: p10, bitwidth: 8)
        markovEntropy = entropy
    }

    // MARK: - File Reading / Queue Fill

    func fillByteQueue(from fileHandle: FileHandle) -> Int {
        var totalLen = 0

        repeat {
            var space = queue.spaceRemaining
            if space > buffSize { space = buffSize }

            guard let data = try? fileHandle.read(upToCount: space), !data.isEmpty else {
                return totalLen
            }

            var buffer = [UInt8](data)
            var len = buffer.count

            // Convert hex to binary if needed
            if opts.hexMode {
                len = hexConverter.convert(&buffer, length: len)
            }

            // Fold uppercase to lowercase
            if opts.fold {
                for i in 0..<len {
                    if buffer[i] >= 0x41 && buffer[i] <= 0x5A {
                        buffer[i] += 32
                    }
                }
            }

            if !opts.wordReverse {
                // Direct transfer to queue
                for i in 0..<len {
                    queue.push(buffer[i], byteReverse: opts.byteReverse)
                    updateMonteCarlo(buffer[i])
                }
                filebytes += UInt64(len)
                totalLen += len
            } else {
                // Word-reverse: accumulate into buffer2, then reverse 4-byte words
                for i in 0..<len {
                    buffer2[buffer2Size] = buffer[i]
                    buffer2Size += 1
                }

                while buffer2Size > 3 {
                    //let base = (buffer2Size > 3) ? buffer2Size - buffer2Size : 0
                    // Process 4 bytes in reverse order
                    for j in 0..<4 {
                        let srcIdx = (totalLen / 4) * 4 + (3 - j)
                        let byte = buffer2[srcIdx < buffer2.count ? srcIdx : 0]
                        let b = opts.byteReverse ? byteReverseTable[Int(byte)] : byte
                        queue.push(b, byteReverse: false)
                        updateMonteCarlo(b)
                    }
                    buffer2Size -= 4
                    totalLen += 4
                }

                if buffer2Size != 0 {
                    // Pad leftover bytes
                    for j in 0..<4 {
                        if j < buffer2Size {
                            queue.push(0x00, byteReverse: false)
                            updateMonteCarlo(0x00)
                        } else {
                            let byte = buffer2[(3 - j)]
                            queue.push(byte, byteReverse: opts.byteReverse)
                            updateMonteCarlo(byte)
                        }
                    }
                    fputs("Warning: Padded \(buffer2Size) extra zeroes for word reverse boundary\n", stderr)
                    buffer2Size = 0
                }
                filebytes += UInt64(len)
            }
        } while queue.spaceRemaining > buffSize

        return totalLen
    }

    // MARK: - Process a single file

    func processFile(_ fileHandle: FileHandle, filename: String, terseIndex: inout Int, parsedFilename: FilenameParsed?) {
        initAll()

        // Skip initial symbols if requested
        if opts.gotSkip {
            if queue.queueUsed == 0 {
                _ = fillByteQueue(from: fileHandle)
            }
            for _ in 0..<opts.skipAmount {
                _ = queue.getSymbol(1)
            }
        }

        // Main analysis loop
        //var running = true
        while true {
            if queue.queueUsed == 0 {
                let bytesRead = fillByteQueue(from: fileHandle)
                if bytesRead == 0 { break }
            }

            let symbol = queue.getSymbol(UInt(opts.symbolLength))
            if symbol == -1 { break }

            symbolCount += 1

            // End if we reach end of substring
            if opts.gotSubstring && (symbolCount + UInt64(opts.skipAmount)) > UInt64(opts.substring) {
                break
            }

            let usym = UInt64(symbol)
            updateMean(usym)
            if !noOccurrenceSpace { updateOccurrences(usym) }
            if !noLongestSpace { updateLongest(usym, symbolPos: symbolCount + UInt64(opts.skipAmount)) }
            updateSCC(usym)
        }

        // Finalize all metrics
        finalizeMean()
        if !noOccurrenceSpace { finalizeOccurrences() }
        if !noLongestSpace { finalizeLongest() }
        finalizeChisq()
        finalizeEntropy()
        finalizeMonteCarlo()
        finalizeCompression()
        finalizeSCC()
        computeMarkov()

        if opts.symbolLength != 8 {
            longestBytePos = (longestPosition * UInt64(opts.symbolLength)) / 8
        } else {
            longestBytePos = longestPosition
        }

        // Output results
        if opts.terse {
            printTerseOutput(terseIndex: terseIndex, filename: filename, parsed: parsedFilename)
        } else {
            printVerboseOutput(filename: filename, parsed: parsedFilename)
        }
    }

    // MARK: - Output

    func printTerseOutput(terseIndex: Int, filename: String, parsed: FilenameParsed?) {
        let r = results

        if opts.entExact {
            let size = opts.symbolLength == 1 ? filebytes * 8 : filebytes
            print("\(terseIndex),\(size),\(r.entropy),\(r.chisqDistribution),\(r.mean),\(r.pi),\(r.scc)")
        } else if let p = parsed, opts.parseFilename, opts.symbolLength == 1 {
            print(String(format: "%4d,%12llu,%8s,%8s,%8.2f,%8.2f,%12f,%12f,%18x,%12f,%12f,%15f,   %16f, %16f, %16f, %16f, %19llx, %18llu, %18f, %15llu, %@",
                         terseIndex, symbolCount, p.deviceID, p.process, p.voltage, p.temperature,
                         r.entropy, r.minEntropy, r.minEntropySymbol, r.chisqPercent, r.mean, r.pi,
                         r.scc, r.p01, r.p10, markovEntropy,
                         longestLongestSymbol, longestLongest, r.longestPValue, longestBytePos, filename))
        } else if !opts.parseFilename && opts.symbolLength == 1 {
            print(String(format: "%4d,%12llu,%11f, %12f,%18x,%12f,%12f,%15f,       %12f, %16f, %16f, %16f, %18llx, %18llu, %18f, %15llu, %@",
                         terseIndex, symbolCount, r.entropy, r.minEntropy, r.minEntropySymbol,
                         r.chisqPercent, r.mean, r.pi, r.scc, r.p01, r.p10, markovEntropy,
                         longestLongestSymbol, longestLongest, r.longestPValue, longestBytePos, filename))
        } else if opts.parseFilename && opts.symbolLength != 1 {
            if let p = parsed {
                print(String(format: "%4d,%12llu,%8s,%8s,%8.2f,%8.2f,%12f,%12f,%18x,%12f,%12f,%15f,   %16f, %16f, %16f, %16f,  %18llx, %18llu,             (null), %15llu, %@",
                             terseIndex, symbolCount, p.deviceID, p.process, p.voltage, p.temperature,
                             r.entropy, r.minEntropy, r.minEntropySymbol, r.chisqPercent, r.mean, r.pi,
                             r.scc, r.p01, r.p10, markovEntropy,
                             longestLongestSymbol, longestLongest, longestBytePos, filename))
            }
        } else {
            print(String(format: "%4d,%12llu,%11f, %12f,%18x,%12f,%12f,%15f,       %12f, %16f, %16f, %16f,  %18llx, %18llu,             (null), %15llu, %@",
                         terseIndex, symbolCount, r.entropy, r.minEntropy, r.minEntropySymbol,
                         r.chisqPercent, r.mean, r.pi, r.scc, r.p01, r.p10, markovEntropy,
                         longestLongestSymbol, longestLongest, longestBytePos, filename))
        }
    }

    func printVerboseOutput(filename: String, parsed: FilenameParsed?) {
        let r = results
        let sl = opts.symbolLength

        // Print occurrence counts if requested
        if opts.printOccurrence && !noOccurrenceSpace {
            for i in 0..<Int(occurrenceSize) {
                let fraction = Double(occurrenceCount[i]) / Double(occurrenceTotal)
                print(String(format: "   Value %4d , frequency=%llu , fraction=%f", i, occurrenceCount[i], fraction))
            }
        }

        // Print longest runs if requested
        if opts.printLongest && !noLongestSpace {
            for i in 0..<Int(occurrenceSize) {
                print(String(format: "   Symbol %x , Longest Run=%llu", i, longestCount[i]))
            }
        }

        if opts.entExact {
            print(String(format: "Entropy = %f bits per byte.\n", r.entropy * (8.0 / Double(sl))))
            print("Optimum compression would reduce the size")
            print(String(format: "of this %d byte file by %d percent\n", Int(filebytes), Int(r.compression)))
            print(String(format: "Chi square distribution for %d samples is %f, and randomly", Int(r.chisqCount), r.chisqDistribution))
            print(String(format: "would exceed this value less than %f percent of the times.\n", r.chisqPercent))
            print(String(format: "Arithmetic mean value of data bytes is %f (127.5 = random).", r.mean))
            print(String(format: "Monte Carlo value for Pi is %f (error %f percent).", r.pi, r.piErr))
            print(String(format: "Serial correlation coefficient is %f (totally uncorrelated = 0.0).", r.scc))
        } else {
            if let p = parsed, opts.parseFilename {
                print("   Device ID   : \(p.deviceID)")
                print("   Process     : \(p.process)")
                print(String(format: "   Voltage     : %0.2fV", p.voltage))
                print(String(format: "   Temperature : %0.2fC", p.temperature))
            }
            print("   Analysing \(symbolCount) \(sl)-bit symbols")
            print(String(format: "   Shannon IID Entropy = %f bits per symbol", r.entropy))
            print(String(format: "   Optimal compression would compress by %f percent", r.compression))
            print(String(format: "   Chi square: symbol count=%llu, distribution=%1.2f, randomly exceeds %1.2f percent of the time",
                         r.chisqCount, r.chisqDistribution, r.chisqPercent))
            print(String(format: "   Mean = %f", r.mean))
            print(String(format: "   Monte Carlo value for Pi is %f (error %1.2f percent).", r.pi, r.piErr))
            print(String(format: "   Serial Correlation = %f", r.scc))
            print(String(format: "   Longest Run Symbol = %llx. Run Length = %llu", longestLongestSymbol, longestLongest))
            if sl == 1 {
                print(String(format: "   Probability of longest run being <= %llu = %f", longestLongest, r.longestPValue))
            }
            print(String(format: "   Position of Longest Run = %llu (0x%llx). Byte position %llu (0x%llx)",
                         longestPosition, longestPosition, longestBytePos, longestBytePos))
            print(String(format: "   A 2 state Markov generator with transition probabilities P01=%f, P10=%f would generate data with entropy %f per bit with 8 bit symbols with the same mean and serial correlation",
                         r.p01, r.p10, markovEntropy))
        }
    }
}

// MARK: - Usage

func displayUsage() {
    let usage = """
    Usage: djent [-brRpcCuhds] [-l <n>] [-i <input file list filename>] [filename] [filename2] ...

Compute statistics of random data.
  Author: David Johnston, dj@deadhat.com

  -i <filename>  --inputfilelist=<filename> Read list of filenames from <filename>
  -p             --parse_filename           Extract CID, Process, Voltage and Temperature from filename.
  -l <n>         --symbol_length=<n>        Treat incoming data symbols as bitlength n. Default is 8.
  -b             --binary                   Treat incoming data as binary. Default bit length will be -l 1
  -r             --byte_reverse             Reverse the bit order in incoming bytes
  -R             --word_reverse             Reverse the byte order in incoming 4 byte words
  -c             --occurrence               Print symbol occurrence counts
  -C             --longest                  Print symbol longest run counts
  -w             --scc_wrap                 Treat data as cyclical in SCC
  -n <n>         --lagn=<n>                 Lag gap in SCC. Default=1
  -S <n>         --skip=<n>                 Skip over <n> initial symbols
  -L <n>         --substring=<n>            Analyse no more than <n> symbols
  -f             --fold                     Fold uppercase letters to lower case
  -t             --terse                    Terse output
  -e             --ent_exact                Exactly match output format of ent
  -s             --suppress_header          Suppress the header in terse output
  -h or -u       --help                     Print this text
"""
  
    let notes = """
Notes
   * By default djent is in hex mode where it reads ascii hex data and converts it to binary to analyze.
     In hex mode, the symbol length defaults to 8, so normal hex files can be treated as a representation
     of bytes. The symbol length can be changed to any value between 1 and 32 bits using the -l <n> option.
   * With the -b option djent switches to binary reads in each byte as binary with a symbol length of 1.
   * To analyze ascii text instead of hex ascii, you need djent to treat each byte as a separate symbol, so
     use binary mode with a symbol length of 8. I.E. djent -b -l 8 <filename>
   * By default djent treats the MSB of each byte as the first. This can be switched so that djent treats
     the LSB as the first bit in each byte using the -r option.
   * Terse output is requested using -t. This outputs in CSV format. The first line is the header. If
     multiple files are provided, there will be one line of CSV output per file in addition to the header.
     The CSV header can be suppressed with -s.
   * To analyze multiple files, just give multiple file names on the command line. To read data in from
     the command line, don't provide a filename and pipe the data in. <datasource> | djent
   * The parse filename option =p picks takes four patterns from the filename to include in the output,
     This is so that it is easy to plot test conditions that are commonly encoded in a filename.
     Fields are delimited by uderscores. The four patters for CID, process, Voltage and Temperature are:
     _CID-<componentID>_ , _PROC-<process info>_, _<x>p<y>V_ and _<x>p<y>C_ . 'p' is the decimal point.
   * To compute the statistics, djent builds a frequency table of the symbols. This can be displayed
     using the -c option. The size of this table is what limits the the maximum symbol size. For each
     of the 2^n symbols, a 64 bit entry in a table is created. So for n=32, that's 32GBytes so the ability
     to handle large symbol sizes is limited by the available memory and the per process allocation limit.
   * The serial correlation coefficient is not wrap around by default, meaning that it does not compare
     the last value in the data with the first. To get wrap around behaviour, use the -w option.
   * The Lag-N correlation coefficient can be computed by using the -n <n> option. This causes the SCC
     computation to compare each Xth symbol with the (X+n)th symbol instead of the (X+1)th symbol.
     If you use wrap around with Lag-N, then the wrap around will reach n bits further into the start
     of the sequence.
   * The byte reverse option -r reverses the order of bits within each byte. The word reverse option -R
     reverses the order of bytes within each 32 bit word, from 3,2,1,0 to 0,1,2,3. Both -R and -r can
     be used together. Using -R with a data that isn't a multiple of 32 bits long will get padded with
     zeros, which may not be what you want. A padding warning will be sent to STDERR.
   * Instead of providing data file names on the command line, djent can be told to read a list of files
     from a text file. The file must have one filename per line. Lines beginning with # will be ignored.
     Use the -i <filename> option to request that djent reads the file list from <filename>.
"""

    let examples = """
Examples
   Print this help
     djent -h

   Analyze hex file from stdin
     cat datafile.hex | djent

   Analyze binary file
     djent -b datafile.bin

   Analyze several files with CSV output
     djent -t data1.hex data2.hex data3.hex

   Analyze ascii symbols - Read in binary and set symbol size to 8.
     djent -b -l 8  textfile.txt

   Analyze binary file with parsable filename.
     djent -b -t -p  rawdata_CID-X23_PROC-TTFT_1p2V_25p0C_.bin
"""

    fputs(usage, stderr)
    print("\n")
    fputs(notes, stderr)
    print("\n")
    fputs(examples, stderr)
    print("\n")
}

// MARK: - Count lines in a file

func countLines(in filename: String) -> Int {
    guard let data = FileManager.default.contents(atPath: filename),
          let text = String(data: data, encoding: .utf8) else {
        return 0
    }
    return text.components(separatedBy: "\n").count
}

// MARK: - Simple getopt-style argument parser

func parseArguments() -> (options: Options, filenames: [String]) {
    var opts = Options()
    var filenames = [String]()
    var gotSymbolLength = false

    let args = CommandLine.arguments
    var i = 1 // skip program name

    while i < args.count {
        let arg = args[i]

        if arg.hasPrefix("-") && arg.count > 1 && !arg.hasPrefix("--") {
            // Short options (can be combined like -bt)
            let chars = Array(arg.dropFirst())
            var ci = 0
            while ci < chars.count {
                let c = chars[ci]
                switch c {
                case "b":
                    if !gotSymbolLength { opts.symbolLength = 1 }
                    opts.hexMode = false
                case "l":
                    i += 1
                    if i < args.count, let v = UInt(args[i]) {
                        opts.symbolLength = v
                        gotSymbolLength = true
                    }
                case "i":
                    i += 1
                    if i < args.count {
                        opts.inputListFilename = args[i]
                        opts.usingInputListFile = true
                    }
                case "c":
                    opts.printOccurrence = true
                case "C":
                    opts.printLongest = true
                case "p":
                    opts.parseFilename = true
                case "r":
                    opts.byteReverse = true
                case "R":
                    opts.wordReverse = true
                case "w":
                    opts.sccWrap = true
                case "n":
                    i += 1
                    if i < args.count, let v = Int(args[i]) { opts.lagN = v }
                case "f":
                    opts.fold = true
                case "t":
                    opts.terse = true
                case "e":
                    opts.entExact = true
                case "s":
                    opts.suppressHeader = true
                case "S":
                    i += 1
                    if i < args.count, let v = Int(args[i]) {
                        opts.gotSkip = true
                        opts.skipAmount = v
                    }
                case "L":
                    i += 1
                    if i < args.count, let v = Int(args[i]) {
                        opts.gotSubstring = true
                        opts.substring = v
                    }
                case "h", "u", "?":
                    displayUsage()
                    exit(0)
                default:
                    break
                }
                ci += 1
            }
        } else if arg.hasPrefix("--") {
            // Long options
            let option = arg.dropFirst(2)
            if option == "help" {
                displayUsage()
                exit(0)
            } else if option == "binary" {
                if !gotSymbolLength { opts.symbolLength = 1 }
                opts.hexMode = false
            } else if option.hasPrefix("symbol_length=") {
                if let v = UInt(option.dropFirst("symbol_length=".count)) {
                    opts.symbolLength = v
                    gotSymbolLength = true
                }
            } else if option.hasPrefix("inputfilelist=") {
                opts.inputListFilename = String(option.dropFirst("inputfilelist=".count))
                opts.usingInputListFile = true
            } else if option == "occurrence" { opts.printOccurrence = true }
            else if option == "longest" { opts.printLongest = true }
            else if option == "parse_filename" { opts.parseFilename = true }
            else if option == "byte_reverse" { opts.byteReverse = true }
            else if option == "word_reverse" { opts.wordReverse = true }
            else if option == "scc_wrap" { opts.sccWrap = true }
            else if option.hasPrefix("lagn=") {
                if let v = Int(option.dropFirst("lagn=".count)) { opts.lagN = v }
            }
            else if option == "fold" { opts.fold = true }
            else if option == "terse" { opts.terse = true }
            else if option == "ent_exact" { opts.entExact = true }
            else if option == "suppress_header" { opts.suppressHeader = true }
            else if option.hasPrefix("skip=") {
                if let v = Int(option.dropFirst("skip=".count)) {
                    opts.gotSkip = true; opts.skipAmount = v
                }
            }
            else if option.hasPrefix("substring=") {
                if let v = Int(option.dropFirst("substring=".count)) {
                    opts.gotSubstring = true; opts.substring = v
                }
            }
        } else {
            // Positional argument = filename
            filenames.append(arg)
        }
        i += 1
    }

    opts.useStdin = filenames.isEmpty && !opts.usingInputListFile
    return (opts, filenames)
}

// MARK: - Main

@main
struct DjentApp {
    static func main() {
        let (opts, filenames) = parseArguments()

        // Validation
        if opts.fold && opts.symbolLength != 8 {
            fputs("Error: Fold must be used with 8 bit word size\n", stderr)
            exit(1)
        }
        if opts.symbolLength < 1 {
            fputs("Error: Symbol length must not be 0 or negative.\n", stderr)
            exit(1)
        }
        if opts.parseFilename && opts.useStdin {
            fputs("Error: Can't parse filename when using stdin for input\n", stderr)
            exit(1)
        }
        if opts.gotSkip && opts.skipAmount < 1 {
            fputs("Error: skip amount must be greater than 0\n", stderr)
            exit(1)
        }
        if opts.gotSubstring && opts.substring < 1 {
            fputs("Error: substring length must be greater than 0\n", stderr)
            exit(1)
        }

        // Build file list
        var fileList: [String]
        if opts.usingInputListFile {
            guard let data = FileManager.default.contents(atPath: opts.inputListFilename),
                  let text = String(data: data, encoding: .utf8) else {
                fputs("Error: Cannot open \(opts.inputListFilename) for reading\n", stderr)
                exit(1)
            }
            fileList = text.components(separatedBy: "\n")
                .map { $0.trimmingCharacters(in: .whitespaces) }
                .filter { !$0.isEmpty && !$0.hasPrefix("#") }
            if fileList.isEmpty {
                fputs("Error: Did not find any filenames in \(opts.inputListFilename)\n", stderr)
                exit(1)
            }
        } else if opts.useStdin {
            fileList = ["<stdin>"]
        } else {
            fileList = filenames
        }

        var terseIndex = 0

        for filename in fileList {
            terseIndex += 1

            let analyzer = EntropyAnalyzer(options: opts)

            // Open file
            let fileHandle: FileHandle
            if filename == "<stdin>" {
                fileHandle = FileHandle.standardInput
            } else {
                // Print the header on first terse file
                if opts.terse && terseIndex == 1 && !opts.suppressHeader {
                    printTerseHeader(opts: opts)
                }

                guard let fh = FileHandle(forReadingAtPath: filename) else {
                    fputs("Error: Unable to open file \(filename)\n", stderr)
                    exit(1)
                }
                fileHandle = fh

                if !opts.terse && !opts.entExact {
                    let mode = opts.hexMode ? "hex text" : "binary"
                    print(" opening \(filename) as \(mode)")
                    print(" Symbol Size(bits) = \(opts.symbolLength)")
                }
            }

            // Print terse header for stdin case
            if filename == "<stdin>" && opts.terse && terseIndex == 1 && !opts.suppressHeader {
                printTerseHeader(opts: opts)
            }

            let parsed: FilenameParsed? = opts.parseFilename ? parseFilename(filename) : nil

            analyzer.processFile(fileHandle, filename: filename, terseIndex: &terseIndex, parsedFilename: parsed)

            if filename != "<stdin>" {
                fileHandle.closeFile()
            }
        }
    }
}

func printTerseHeader(opts: Options) {
    if opts.entExact {
        if opts.symbolLength == 1 {
            print("0,File-bits,Entropy,Chi-square,Mean,Monte-Carlo-Pi,Serial-Correlation")
        } else {
            print("0,File-bytes,Entropy,Chi-square,Mean,Monte-Carlo-Pi,Serial-Correlation")
        }
    } else if opts.parseFilename {
        print("   0,     symbols,     CID, Process, Voltage,    Temp,     Entropy,  MinEntropy, MinEntropy-Symbol,  Chi-square,        Mean, Monte-Carlo-Pi, Serial-Correlation,              P01,              P10,  mkv_min_entropy,  Longest-Run-Symbol, Longest-Run-Length, Longest-Run-PValue, Longest-Run-Pos, Filename")
    } else {
        print("   0,     symbols,    Entropy,  Min_entropy, MinEntropy-Symbol,  Chi-square,        Mean, Monte-Carlo-Pi, Serial-Correlation,              P01,              P10,  mkv_min_entropy, Longest-Run-Symbol, Longest-Run-Length, Longest-Run-PValue, Longest-Run-Pos, Filename")
    }
}
