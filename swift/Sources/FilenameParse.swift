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

/// Results from parsing a filename for test conditions
struct FilenameParsed {
    var voltage: Double = 0.0
    var temperature: Double = 0.0
    var deviceID: String = ""
    var process: String = ""
}

/// Look for the _<int>p<int>V_ voltage pattern in a string.
/// Returns the matched substring or nil.
func findVoltagePattern(in str: String) -> String? {
    let chars = Array(str)
    let len = chars.count
    var start = 0
    var state = 1
    var pos = 0

    while pos < len {
        let c = chars[pos]
        switch state {
        case 1: // looking for _
            if c == "_" { state = 2; start = pos }
            pos += 1
        case 2: // first digit
            if c.isNumber { state = 3 } else { state = 1 }
            pos += 1
        case 3: // rest of first int
            if c.isNumber { /* stay */ }
            else if c == "p" || c == "." { state = 4 }
            else { state = 1 }
            pos += 1
        case 4: // first digit after decimal
            if c.isNumber { state = 5 } else { state = 1 }
            pos += 1
        case 5: // rest of second int
            if c.isNumber { /* stay */ }
            else if c == "V" { state = 6 }
            else { state = 1 }
            pos += 1
        case 6: // trailing _ or .
            if c == "_" || c == "." {
                let end = pos
                return String(chars[start...end])
            } else { state = 1 }
            pos += 1
        default:
            pos += 1
        }
    }
    return nil
}

/// Look for the _<int>p<int>C_ temperature pattern in a string.
func findTemperaturePattern(in str: String) -> String? {
    let chars = Array(str)
    let len = chars.count
    var start = 0
    var state = 1
    var pos = 0

    while pos < len {
        let c = chars[pos]
        switch state {
        case 1:
            if c == "_" { state = 2; start = pos }
            pos += 1
        case 2: // first digit or minus sign
            if c.isNumber || c == "-" { state = 3 } else { state = 1 }
            pos += 1
        case 3: // rest of first int
            if c.isNumber { /* stay */ }
            else if c == "p" || c == "." { state = 4 }
            else { state = 1 }
            pos += 1
        case 4:
            if c.isNumber { state = 5 } else { state = 1 }
            pos += 1
        case 5:
            if c.isNumber { /* stay */ }
            else if c == "C" { state = 6 }
            else { state = 1 }
            pos += 1
        case 6:
            if c == "_" || c == "." {
                let end = pos
                return String(chars[start...end])
            } else { state = 1 }
            pos += 1
        default:
            pos += 1
        }
    }
    return nil
}

/// Look for the _CID-<ID>_ pattern in a string.
func findCIDPattern(in str: String) -> String? {
    let chars = Array(str)
    let len = chars.count
    var start = 0
    var state = 1
    var pos = 0

    while pos < len {
        let c = chars[pos]
        switch state {
        case 1:
            if c == "_" { state = 2; start = pos }
            pos += 1
        case 2:
            state = c == "C" ? 3 : 1; pos += 1
        case 3:
            state = c == "I" ? 4 : 1; pos += 1
        case 4:
            state = c == "D" ? 5 : 1; pos += 1
        case 5:
            state = c == "-" ? 6 : 1; pos += 1
        case 6: // first char of ID
            if c != "_" { state = 7 } else { state = 1 }
            pos += 1
        case 7: // rest of ID
            if c == "_" {
                let end = pos
                return String(chars[start...end])
            }
            pos += 1
        default:
            pos += 1
        }
    }
    return nil
}

/// Look for the _PROC-<name>_ pattern in a string.
func findProcPattern(in str: String) -> String? {
    let chars = Array(str)
    let len = chars.count
    var start = 0
    var state = 1
    var pos = 0

    while pos < len {
        let c = chars[pos]
        switch state {
        case 1:
            if c == "_" { state = 2; start = pos }
            pos += 1
        case 2:
            state = c == "P" ? 3 : 1; pos += 1
        case 3:
            state = c == "R" ? 4 : 1; pos += 1
        case 4:
            state = c == "O" ? 5 : 1; pos += 1
        case 5:
            state = c == "C" ? 6 : 1; pos += 1
        case 6:
            state = c == "-" ? 7 : 1; pos += 1
        case 7:
            if c != "_" { state = 8 } else { state = 1 }
            pos += 1
        case 8:
            if c == "_" {
                let end = pos
                return String(chars[start...end])
            }
            pos += 1
        default:
            pos += 1
        }
    }
    return nil
}

/// Parse a filename for voltage, temperature, CID and process fields.
func parseFilename(_ filename: String) -> FilenameParsed {
    var result = FilenameParsed()

    if var match = findVoltagePattern(in: filename) {
        match = match.replacingOccurrences(of: "p", with: ".")
        // Remove leading _ and trailing V_ or V.
        let trimmed = match.trimmingCharacters(in: CharacterSet(charactersIn: "_."))
        if let vStr = trimmed.components(separatedBy: "V").first,
           let v = Double(vStr) {
            result.voltage = v
        }
    } else {
        fputs("Regex error scanning for _<num>p<num>V_:\n", stderr)
    }

    if var match = findTemperaturePattern(in: filename) {
        match = match.replacingOccurrences(of: "p", with: ".")
        let trimmed = match.trimmingCharacters(in: CharacterSet(charactersIn: "_."))
        if let tStr = trimmed.components(separatedBy: "C").first,
           let t = Double(tStr) {
            result.temperature = t
        }
    } else {
        fputs("Regex error scanning for _<num>p<num>C_:\n", stderr)
    }

    if let match = findCIDPattern(in: filename) {
        // Extract ID after _CID- and before trailing _
        let inner = match.dropFirst(5).dropLast(1) // drop "_CID-" and trailing "_"
        result.deviceID = String(inner)
    } else {
        fputs("Regex error scanning for _CID-<ID>_:\n", stderr)
    }

    if let match = findProcPattern(in: filename) {
        let inner = match.dropFirst(6).dropLast(1) // drop "_PROC-" and trailing "_"
        result.process = String(inner)
    } else {
        fputs("Regex error scanning for _PROC-<ID>_:\n", stderr)
    }

    return result
}
