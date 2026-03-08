//! djent - Filename parsing for test condition extraction.
//!
//! Copyright (C) 2017  David Johnston
//! Translated to Zig 0.15.x, 2025
//!
//! Contact. David Johnston dj@deadhat.com

const std = @import("std");

pub const FilenameParsed = struct {
    voltage: f64 = 0.0,
    temperature: f64 = 0.0,
    device_id: [256]u8 = [_]u8{0} ** 256,
    device_id_len: usize = 0,
    process: [256]u8 = [_]u8{0} ** 256,
    process_len: usize = 0,
};

const MatchResult = struct { start: usize, end: usize };

/// Look for _<int>p<int>V_ voltage pattern. Returns matched range or null.
pub fn findVPattern(str: []const u8) ?MatchResult {
    var state: u8 = 1;
    var pos: usize = 0;
    var start: usize = 0;

    while (pos < str.len) {
        const c = str[pos];
        switch (state) {
            1 => {
                if (c == '_') { state = 2; start = pos; }
                pos += 1;
            },
            2 => {
                state = if (std.ascii.isDigit(c)) 3 else 1;
                pos += 1;
            },
            3 => {
                if (std.ascii.isDigit(c)) {
                    // stay
                } else if (c == 'p' or c == '.') {
                    state = 4;
                } else { state = 1; }
                pos += 1;
            },
            4 => {
                state = if (std.ascii.isDigit(c)) 5 else 1;
                pos += 1;
            },
            5 => {
                if (std.ascii.isDigit(c)) {
                    // stay
                } else if (c == 'V') { state = 6; } else { state = 1; }
                pos += 1;
            },
            6 => {
                if (c == '_' or c == '.') return .{ .start = start, .end = pos };
                state = 1;
                pos += 1;
            },
            else => { pos += 1; },
        }
    }
    return null;
}

/// Look for _<int>p<int>C_ temperature pattern.
pub fn findTPattern(str: []const u8) ?MatchResult {
    var state: u8 = 1;
    var pos: usize = 0;
    var start: usize = 0;

    while (pos < str.len) {
        const c = str[pos];
        switch (state) {
            1 => {
                if (c == '_') { state = 2; start = pos; }
                pos += 1;
            },
            2 => {
                state = if (std.ascii.isDigit(c) or c == '-') 3 else 1;
                pos += 1;
            },
            3 => {
                if (std.ascii.isDigit(c)) {} else if (c == 'p' or c == '.') { state = 4; } else { state = 1; }
                pos += 1;
            },
            4 => {
                state = if (std.ascii.isDigit(c)) 5 else 1;
                pos += 1;
            },
            5 => {
                if (std.ascii.isDigit(c)) {} else if (c == 'C') { state = 6; } else { state = 1; }
                pos += 1;
            },
            6 => {
                if (c == '_' or c == '.') return .{ .start = start, .end = pos };
                state = 1;
                pos += 1;
            },
            else => { pos += 1; },
        }
    }
    return null;
}

/// Look for _CID-<ID>_ pattern.
pub fn findCIDPattern(str: []const u8) ?MatchResult {
    var state: u8 = 1;
    var pos: usize = 0;
    var start: usize = 0;

    while (pos < str.len) {
        const c = str[pos];
        switch (state) {
            1 => { if (c == '_') { state = 2; start = pos; } pos += 1; },
            2 => { state = if (c == 'C') 3 else 1; pos += 1; },
            3 => { state = if (c == 'I') 4 else 1; pos += 1; },
            4 => { state = if (c == 'D') 5 else 1; pos += 1; },
            5 => { state = if (c == '-') 6 else 1; pos += 1; },
            6 => { state = if (c != '_') 7 else 1; pos += 1; },
            7 => {
                if (c == '_') return .{ .start = start, .end = pos };
                pos += 1;
            },
            else => { pos += 1; },
        }
    }
    return null;
}

/// Look for _PROC-<n>_ pattern.
pub fn findProcPattern(str: []const u8) ?MatchResult {
    var state: u8 = 1;
    var pos: usize = 0;
    var start: usize = 0;

    while (pos < str.len) {
        const c = str[pos];
        switch (state) {
            1 => { if (c == '_') { state = 2; start = pos; } pos += 1; },
            2 => { state = if (c == 'P') 3 else 1; pos += 1; },
            3 => { state = if (c == 'R') 4 else 1; pos += 1; },
            4 => { state = if (c == 'O') 5 else 1; pos += 1; },
            5 => { state = if (c == 'C') 6 else 1; pos += 1; },
            6 => { state = if (c == '-') 7 else 1; pos += 1; },
            7 => { state = if (c != '_') 8 else 1; pos += 1; },
            8 => {
                if (c == '_') return .{ .start = start, .end = pos };
                pos += 1;
            },
            else => { pos += 1; },
        }
    }
    return null;
}

/// Parse a float from a matched pattern like "_1p2V_" -> 1.2
fn parseMatchedFloat(str: []const u8, suffix: u8) f64 {
    var buf: [256]u8 = undefined;
    var len: usize = 0;
    var started = false;

    for (str) |c| {
        if (c == '_' or c == '.') {
            if (started) break;
            continue;
        }
        if (c == suffix) break;
        started = true;
        buf[len] = if (c == 'p') '.' else c;
        len += 1;
    }

    return std.fmt.parseFloat(f64, buf[0..len]) catch 0.0;
}

/// Parse a filename for voltage, temperature, CID and process fields.
pub fn parseFilename(filename: []const u8) FilenameParsed {
    var result = FilenameParsed{};

    if (findVPattern(filename)) |m| {
        result.voltage = parseMatchedFloat(filename[m.start .. m.end + 1], 'V');
    } else {
        std.debug.print("Regex error scanning for _<num>p<num>V_:\n", .{});
    }

    if (findTPattern(filename)) |m| {
        result.temperature = parseMatchedFloat(filename[m.start .. m.end + 1], 'C');
    } else {
        std.debug.print("Regex error scanning for _<num>p<num>C_:\n", .{});
    }

    if (findCIDPattern(filename)) |m| {
        const id_start = m.start + 5; // skip "_CID-"
        const id_end = m.end;
        const id_len = id_end - id_start;
        if (id_len > 0 and id_len <= 255) {
            @memcpy(result.device_id[0..id_len], filename[id_start..id_end]);
            result.device_id_len = id_len;
        }
    } else {
        std.debug.print("Regex error scanning for _CID-<ID>_:\n", .{});
    }

    if (findProcPattern(filename)) |m| {
        const id_start = m.start + 6; // skip "_PROC-"
        const id_end = m.end;
        const id_len = id_end - id_start;
        if (id_len > 0 and id_len <= 255) {
            @memcpy(result.process[0..id_len], filename[id_start..id_end]);
            result.process_len = id_len;
        }
    } else {
        std.debug.print("Regex error scanning for _PROC-<ID>_:\n", .{});
    }

    return result;
}
