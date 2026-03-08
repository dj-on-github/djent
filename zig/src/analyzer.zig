//! djent - Core entropy analysis engine.
//!
//! Copyright (C) 2017  David Johnston
//! Translated to Zig 0.15.x, 2025
//!
//! Contact. David Johnston dj@deadhat.com

const std = @import("std");
const math = std.math;
const mathy = @import("mathy_things.zig");
const markov = @import("markov2p.zig");
const longest_run = @import("longest_run_cdf.zig");
const fnparse = @import("filename_parse.zig");

const queue_capacity = 4096;
const buff_size = 2048;

// Simple stdout helper for Zig 0.15 (bufPrint + writeAll).
fn stdoutPrint(comptime fmt: []const u8, args: anytype) void {
    var buf: [8192]u8 = undefined;
    const msg = std.fmt.bufPrint(&buf, fmt, args) catch return;
    const f = std.fs.File.stdout();
    f.writeAll(msg) catch {};
}

pub const byte_reverse_table = [256]u8{
    0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0,
    0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8,
    0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4,
    0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC,
    0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2,
    0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
    0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6,
    0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
    0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
    0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9,
    0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
    0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
    0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3,
    0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
    0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7,
    0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF,
};

pub const Options = struct {
    symbol_length: u32 = 8,
    hex_mode: bool = true,
    print_occurrence: bool = false,
    print_longest: bool = false,
    fold: bool = false,
    terse: bool = false,
    use_stdin: bool = true,
    scc_wrap: bool = false,
    lag_n: u32 = 1,
    using_input_list_file: bool = false,
    input_list_filename: []const u8 = "",
    suppress_header: bool = false,
    byte_reverse: bool = false,
    parse_filename: bool = false,
    word_reverse: bool = false,
    ent_exact: bool = false,
    got_skip: bool = false,
    skip_amount: u64 = 0,
    got_substring: bool = false,
    substring: u64 = 0,
};

pub const Results = struct {
    mean: f64 = 0.0,
    chisq_count: u64 = 0,
    chisq_distribution: f64 = 0.0,
    chisq_percent: f64 = 0.0,
    entropy: f64 = 0.0,
    min_entropy: f64 = 0.0,
    min_entropy_symbol: u32 = 0,
    pi: f64 = 0.0,
    pi_err: f64 = 0.0,
    compression: f64 = 0.0,
    scc: f64 = 0.0,
    p01: f64 = 0.0,
    p10: f64 = 0.0,
    longest_pvalue: f64 = 0.0,
};

// MARK: - Byte Queue

pub const ByteQueue = struct {
    queue: [queue_capacity]u8 = [_]u8{0} ** queue_capacity,
    q_start: usize = 0,
    q_end: usize = 0,
    q_used: usize = 0,
    current_byte: u32 = 0,
    bits_used: u32 = 0,
    got_byte: bool = false,
    current_symbol: i64 = 0,
    bits_in_symbol: u32 = 0,
    symbol_mask: u64 = 0,

    pub fn init(self: *ByteQueue, symbol_length: u32) void {
        self.q_start = 0;
        self.q_end = 0;
        self.q_used = 0;
        self.got_byte = false;
        self.current_byte = 0;
        self.bits_used = 0;
        self.current_symbol = 0;
        self.bits_in_symbol = 0;
        @memset(&self.queue, 0);
        self.symbol_mask = mathy.ipow(2, symbol_length) -% 1;
    }

    pub fn push(self: *ByteQueue, byte: u8, do_reverse: bool) void {
        const b = if (do_reverse) byte_reverse_table[byte] else byte;
        self.queue[self.q_end] = b;
        self.q_end = (self.q_end + 1) % queue_capacity;
        self.q_used += 1;
    }

    pub fn spaceRemaining(self: *const ByteQueue) usize {
        return queue_capacity - self.q_used;
    }

    pub fn getSymbol(self: *ByteQueue, symbol_length: u32) i64 {
        self.current_symbol = 0;
        if (!self.got_byte) {
            if (self.q_used * 8 < symbol_length) return -1;
            self.current_byte = self.queue[self.q_start];
            self.q_start = (self.q_start + 1) % queue_capacity;
            self.q_used -= 1;
            self.bits_used = 0;
            self.got_byte = true;
        }
        if (symbol_length == 1) {
            self.current_symbol = @intCast((self.current_byte & 0x80) >> 7);
            self.current_byte = (self.current_byte << 1) & 0xFF;
            self.bits_used += 1;
            if (self.bits_used == 8) { self.got_byte = false; self.bits_used = 0; }
            return self.current_symbol;
        }
        if (symbol_length == 8) {
            self.current_symbol = @intCast(self.current_byte);
            self.got_byte = false;
            return self.current_symbol;
        }
        self.bits_in_symbol = 0;
        while (true) {
            const temp: u64 = (self.current_byte & 0x80) >> 7;
            self.current_byte = (self.current_byte << 1) & 0xFF;
            self.bits_used += 1;
            if (self.bits_used == 8) { self.got_byte = false; self.bits_used = 0; }
            const shifted: u64 = (@as(u64, @bitCast(self.current_symbol)) << 1) | temp;
            self.current_symbol = @bitCast(shifted & self.symbol_mask);
            self.bits_in_symbol += 1;
            if (!self.got_byte) {
                if (self.q_used * 8 < symbol_length) return -1;
                self.current_byte = self.queue[self.q_start];
                self.q_start = (self.q_start + 1) % queue_capacity;
                self.q_used -= 1;
                self.bits_used = 0;
                self.got_byte = true;
            }
            if (self.bits_in_symbol >= symbol_length) break;
        }
        return self.current_symbol;
    }
};

// MARK: - Hex Converter

pub const HexConverter = struct {
    hexpair: [2]u8 = .{ 0, 0 },
    state: u8 = 0,

    pub fn reset(self: *HexConverter) void {
        self.hexpair = .{ 0, 0 };
        self.state = 0;
    }

    fn isHex(c: u8) bool {
        return (c >= '0' and c <= '9') or (c >= 'A' and c <= 'F') or (c >= 'a' and c <= 'f');
    }

    fn isHexOrX(c: u8) bool {
        return isHex(c) or c == 'x';
    }

    pub fn convert(self: *HexConverter, buffer: []u8, length: usize) usize {
        var out_pos: usize = 0;
        var scan_pos: usize = 0;
        while (scan_pos < length) {
            if (self.state == 0) {
                if (isHex(buffer[scan_pos])) { self.hexpair[0] = buffer[scan_pos]; self.state = 1; }
                scan_pos += 1;
            } else if (self.state == 1) {
                const c = buffer[scan_pos];
                if ((isHexOrX(c) and self.hexpair[0] == '0') or isHex(c)) { self.hexpair[1] = c; self.state = 2; }
                scan_pos += 1;
            } else if (self.state == 2) {
                if (self.hexpair[0] == '0' and self.hexpair[1] == 'x') {
                    self.state = 0;
                } else {
                    var nybble: u8 = 0;
                    const h0 = self.hexpair[0];
                    if (h0 >= '0' and h0 <= '9') { nybble = h0 - '0'; } else if (h0 >= 'A' and h0 <= 'F') { nybble = h0 - 'A' + 10; } else if (h0 >= 'a' and h0 <= 'f') { nybble = h0 - 'a' + 10; }
                    nybble <<= 4;
                    var byte_val: u8 = nybble;
                    const h1 = self.hexpair[1];
                    if (h1 >= '0' and h1 <= '9') { byte_val += h1 - '0'; } else if (h1 >= 'A' and h1 <= 'F') { byte_val += h1 - 'A' + 10; } else if (h1 >= 'a' and h1 <= 'f') { byte_val += h1 - 'a' + 10; }
                    buffer[out_pos] = byte_val;
                    out_pos += 1;
                    self.state = 0;
                }
            }
        }
        return out_pos;
    }
};

// MARK: - Analyzer

pub const Analyzer = struct {
    const Self = @This();

    opts: Options,
    queue: ByteQueue = .{},
    hex_conv: HexConverter = .{},
    symbol_count: u64 = 0,
    mean_total: u64 = 0,
    filebytes: u64 = 0,
    ent: f64 = 0.0,
    occurrence_size: u64 = 0,
    occurrence_count: []u64 = &.{},
    occurrence_total: u64 = 0,
    longest_count: []u64 = &.{},
    longest_last_symbol: u64 = 0,
    longest_run: u64 = 0,
    longest_longest: u64 = 0,
    longest_longest_symbol: u64 = 0,
    longest_position: u64 = 0,
    longest_new_pos: u64 = 0,
    longest_byte_pos: u64 = 0,
    chisq: f64 = 0.0,
    chisq_sum: f64 = 0.0,
    chisq_prob: []f64 = &.{},
    mp: u64 = 0,
    monty_total: u64 = 0,
    monty_inside: u64 = 0,
    radius_sq: f64 = 0.0,
    position_x: f64 = 0.0,
    position_y: f64 = 0.0,
    monte: [6]u64 = .{ 0, 0, 0, 0, 0, 0 },
    t1: u64 = 0, t2: u64 = 0, t3: u64 = 0,
    scc_first: bool = true,
    first_symbol: u64 = 0,
    scc_previous: u64 = 0,
    scc_count: u64 = 0,
    scc_fifo: [256]u64 = [_]u64{0} ** 256,
    scc_first_lagn: [256]u64 = [_]u64{0} ** 256,
    aeqb_count: u64 = 0,
    mean_count: u64 = 0,
    results: Results = .{},
    markov_entropy: f64 = 0.0,
    buffer2: [buff_size + 4]u8 = [_]u8{0} ** (buff_size + 4),
    buffer2_size: usize = 0,
    allocator: std.mem.Allocator,

    pub fn create(allocator: std.mem.Allocator, opts: Options) Self {
        return Self{ .opts = opts, .allocator = allocator };
    }

    pub fn initAll(self: *Self) !void {
        self.symbol_count = 0;
        self.filebytes = 0;
        self.buffer2_size = 0;
        self.queue.init(self.opts.symbol_length);
        self.hex_conv.reset();
        self.mean_total = 0;
        self.ent = 0.0;
        self.occurrence_total = 0;
        self.occurrence_size = mathy.ipow(2, self.opts.symbol_length);
        if (self.occurrence_count.len > 0) self.allocator.free(self.occurrence_count);
        self.occurrence_count = try self.allocator.alloc(u64, @intCast(self.occurrence_size));
        @memset(self.occurrence_count, 0);
        if (self.longest_count.len > 0) self.allocator.free(self.longest_count);
        self.longest_count = try self.allocator.alloc(u64, @intCast(self.occurrence_size));
        @memset(self.longest_count, 0);
        self.longest_last_symbol = 0; self.longest_run = 0; self.longest_longest = 0;
        self.longest_longest_symbol = 0; self.longest_position = 0; self.longest_new_pos = 0; self.longest_byte_pos = 0;
        self.chisq = 0.0; self.chisq_sum = 0.0;
        if (self.chisq_prob.len > 0) self.allocator.free(self.chisq_prob);
        self.chisq_prob = try self.allocator.alloc(f64, @intCast(self.occurrence_size));
        @memset(self.chisq_prob, 0.0);
        self.mp = 0; self.monty_total = 0; self.monty_inside = 0;
        self.radius_sq = (256.0 * 256.0 * 256.0) - 1.0; self.radius_sq *= self.radius_sq;
        self.monte = .{ 0, 0, 0, 0, 0, 0 };
        self.t1 = 0; self.t2 = 0; self.t3 = 0;
        self.scc_first = true; self.scc_previous = 0; self.scc_count = 0; self.first_symbol = 0;
        @memset(&self.scc_fifo, 0); @memset(&self.scc_first_lagn, 0);
        self.aeqb_count = 0; self.mean_count = 0;
        self.results = .{};
        self.markov_entropy = 0.0;
    }

    pub fn deinit(self: *Self) void {
        if (self.occurrence_count.len > 0) self.allocator.free(self.occurrence_count);
        if (self.longest_count.len > 0) self.allocator.free(self.longest_count);
        if (self.chisq_prob.len > 0) self.allocator.free(self.chisq_prob);
        self.occurrence_count = &.{};
        self.longest_count = &.{};
        self.chisq_prob = &.{};
    }

    // -- Update routines --
    pub fn updateMean(self: *Self, symbol: u64) void { self.mean_total +%= symbol; }
    pub fn updateOccurrences(self: *Self, symbol: u64) void { self.occurrence_count[@intCast(symbol)] += 1; self.occurrence_total += 1; }

    pub fn updateLongest(self: *Self, symbol: u64, symbol_pos: u64) void {
        const symbol_byte_pos = (@as(u64, self.opts.symbol_length) * symbol_pos) / 8;
        if (symbol == self.longest_last_symbol) {
            self.longest_run += 1;
            if (self.longest_run > self.longest_count[@intCast(symbol)]) self.longest_count[@intCast(symbol)] = self.longest_run;
            if (self.longest_run > self.longest_longest) { self.longest_longest = self.longest_run; self.longest_longest_symbol = symbol; self.longest_position = self.longest_new_pos; }
        } else { self.longest_run = 1; self.longest_last_symbol = symbol; self.longest_new_pos = symbol_byte_pos; }
    }

    pub fn updateMonteCarlo(self: *Self, symbol: u8) void {
        self.monte[@intCast(self.mp)] = symbol;
        self.mp += 1;
        if (self.mp > 5) {
            self.mp = 0; self.monty_total += 1;
            self.position_x = 0.0; self.position_y = 0.0;
            for (0..3) |mj| {
                self.position_x = (self.position_x * 256.0) + @as(f64, @floatFromInt(self.monte[mj]));
                self.position_y = (self.position_y * 256.0) + @as(f64, @floatFromInt(self.monte[3 + mj]));
            }
            if ((self.position_x * self.position_x) + (self.position_y * self.position_y) <= self.radius_sq) self.monty_inside += 1;
        }
    }

    pub fn updateSCC(self: *Self, symbol: u64) void {
        const lagn = self.opts.lag_n;
        if (lagn == 1) {
            self.scc_count += 1;
            if (self.scc_first) { self.scc_first = false; self.first_symbol = symbol; } else { self.t1 +%= self.scc_previous *% symbol; if (self.scc_previous == symbol) self.aeqb_count += 1; }
            self.mean_count +%= symbol; self.t2 +%= symbol *% symbol; self.t3 +%= symbol; self.scc_previous = symbol;
        } else {
            self.scc_count += 1;
            if (self.scc_count <= lagn) { self.scc_fifo[@intCast(self.scc_count - 1)] = symbol; self.scc_first_lagn[@intCast(self.scc_count - 1)] = symbol; } else {
                self.t1 +%= self.scc_fifo[0] *% symbol; if (self.scc_fifo[0] == symbol) self.aeqb_count += 1;
                for (0..lagn) |i| { self.scc_fifo[i] = self.scc_fifo[i + 1]; }
                self.mean_count +%= symbol; self.scc_fifo[lagn] = symbol; self.t2 +%= symbol *% symbol; self.t3 +%= symbol;
            }
        }
    }

    // -- Finalization routines --
    pub fn finalizeMean(self: *Self) void { self.results.mean = @as(f64, @floatFromInt(self.mean_total)) / @as(f64, @floatFromInt(self.symbol_count)); }

    pub fn finalizeEntropy(self: *Self) void {
        self.ent = 0.0;
        for (0..@intCast(self.occurrence_size)) |i| { if (self.chisq_prob[i] > 0.0) self.ent += self.chisq_prob[i] * @log10(1.0 / self.chisq_prob[i]) * 3.32192809488736234787; }
        self.results.entropy = self.ent;
    }

    pub fn finalizeOccurrences(self: *Self) void {
        var max_c: u64 = 0; var max_symbol: u32 = 0;
        for (0..@intCast(self.occurrence_size)) |i| { if (self.occurrence_count[i] > max_c) { max_c = self.occurrence_count[i]; max_symbol = @intCast(i); } }
        const max_p = @as(f64, @floatFromInt(max_c)) / @as(f64, @floatFromInt(self.occurrence_total));
        const max_p_ent = (-@log10(max_p) / @log10(@as(f64, 2.0))) / @as(f64, @floatFromInt(self.opts.symbol_length));
        self.results.min_entropy = max_p_ent;
        self.results.min_entropy_symbol = max_symbol;
        if (!self.opts.terse) stdoutPrint("   Min Entropy (by max occurrence of symbol {x}) = {d:.6}\n", .{ max_symbol, max_p_ent });
    }

    pub fn finalizeLongest(self: *Self) void {
        self.results.longest_pvalue = longest_run.longestRunCDF(@intCast(self.longest_longest), self.symbol_count);
        if (self.opts.symbol_length != 8) { self.longest_byte_pos = (self.symbol_count * self.opts.symbol_length) / 8; } else { self.longest_byte_pos = self.symbol_count; }
    }

    pub fn finalizeChisq(self: *Self) void {
        const expected = @as(f64, @floatFromInt(self.occurrence_total)) / @as(f64, @floatFromInt(self.occurrence_size));
        self.chisq = 0.0; self.chisq_sum = 0.0;
        for (0..@intCast(self.occurrence_size)) |i| {
            const diff = @as(f64, @floatFromInt(self.occurrence_count[i])) - expected;
            self.chisq_prob[i] = @as(f64, @floatFromInt(self.occurrence_count[i])) / @as(f64, @floatFromInt(self.occurrence_total));
            self.chisq += (diff * diff) / expected;
            self.chisq_sum += @as(f64, @floatFromInt(i * self.occurrence_count[i]));
        }
        self.results.chisq_count = self.occurrence_total;
        self.results.chisq_distribution = self.chisq;
        self.results.chisq_percent = mathy.chisqp(self.chisq, @intCast(self.occurrence_size - 1)) * 100.0;
    }

    pub fn finalizeMonteCarlo(self: *Self) void {
        if (self.monty_total == 0) { self.results.pi = 0.0; self.results.pi_err = 100.0; return; }
        const montepi = 4.0 * (@as(f64, @floatFromInt(self.monty_inside)) / @as(f64, @floatFromInt(self.monty_total)));
        self.results.pi = montepi;
        self.results.pi_err = (@abs(math.pi - montepi) / math.pi) * 100.0;
    }

    pub fn finalizeCompression(self: *Self) void {
        const sl: f64 = @floatFromInt(self.opts.symbol_length);
        self.results.compression = (100.0 * (sl - self.ent)) / sl;
    }

    pub fn finalizeSCC(self: *Self) void {
        const lagn = self.opts.lag_n;
        if (self.opts.scc_wrap) {
            if (lagn == 1) { self.t1 +%= self.scc_previous *% self.first_symbol; self.t2 +%= self.first_symbol *% self.first_symbol; self.t3 +%= self.first_symbol; } else { for (0..lagn) |i| { self.t1 +%= self.scc_fifo[i] *% self.scc_first_lagn[i]; self.t2 +%= self.scc_first_lagn[i] *% self.scc_first_lagn[i]; self.t3 +%= self.scc_first_lagn[i]; } }
        } else { self.scc_count -= lagn; }
        const top: i64 = @bitCast(self.scc_count *% self.t1 -% self.t3 *% self.t3);
        const bottom: i64 = @bitCast(self.scc_count *% self.t2 -% self.t3 *% self.t3);
        self.results.scc = @as(f64, @floatFromInt(top)) / @as(f64, @floatFromInt(bottom));
    }

    pub fn computeMarkov(self: *Self) void {
        const p01 = self.results.mean * (1.0 - self.results.scc);
        const p10 = (1.0 - self.results.mean) * (1.0 - self.results.scc);
        self.results.p01 = p01; self.results.p10 = p10;
        self.markov_entropy = markov.pToEntropy(p01, p10, 8).entropy;
    }

    // -- File reading (uses std.fs.File directly) --
    pub fn fillByteQueue(self: *Self, file: std.fs.File) !usize {
        var total_len: usize = 0;
        while (true) {
            var space = self.queue.spaceRemaining();
            if (space > buff_size) space = buff_size;
            if (space == 0) break;
            var buffer: [buff_size]u8 = undefined;
            const bytes_read = try file.read(buffer[0..space]);
            if (bytes_read == 0) return total_len;
            var len = bytes_read;
            if (self.opts.hex_mode) len = self.hex_conv.convert(&buffer, len);
            if (self.opts.fold) { for (0..len) |i| { buffer[i] = std.ascii.toLower(buffer[i]); } }
            if (!self.opts.word_reverse) {
                for (0..len) |i| { self.queue.push(buffer[i], self.opts.byte_reverse); self.updateMonteCarlo(buffer[i]); }
                self.filebytes += len; total_len += len;
            } else {
                for (0..len) |i| { self.buffer2[self.buffer2_size] = buffer[i]; self.buffer2_size += 1; }
                var wi: usize = 0;
                while (self.buffer2_size > 3) {
                    for (0..4) |j| { const src_idx = wi * 4 + (3 - j); const byte = self.buffer2[src_idx]; if (self.opts.byte_reverse) { self.queue.push(byte_reverse_table[byte], false); self.updateMonteCarlo(byte_reverse_table[byte]); } else { self.queue.push(byte, false); self.updateMonteCarlo(byte); } }
                    self.buffer2_size -= 4; total_len += 4; wi += 1;
                }
                if (self.buffer2_size != 0) { std.debug.print("Warning: Padded {d} extra zeroes for word reverse boundary\n", .{self.buffer2_size}); self.buffer2_size = 0; }
                self.filebytes += len;
            }
            if (self.queue.spaceRemaining() <= buff_size) break;
        }
        return total_len;
    }

    /// Process data from a file handle.
    pub fn process(self: *Self, file: std.fs.File) !void {
        try self.initAll();
        if (self.opts.got_skip) {
            if (self.queue.q_used == 0) _ = try self.fillByteQueue(file);
            for (0..@intCast(self.opts.skip_amount)) |_| { _ = self.queue.getSymbol(1); }
        }
        while (true) {
            if (self.queue.q_used == 0) { if (try self.fillByteQueue(file) == 0) break; }
            const symbol = self.queue.getSymbol(self.opts.symbol_length);
            if (symbol == -1) break;
            self.symbol_count += 1;
            if (self.opts.got_substring and (self.symbol_count + self.opts.skip_amount) > self.opts.substring) break;
            const usym: u64 = @intCast(symbol);
            self.updateMean(usym);
            self.updateOccurrences(usym);
            self.updateLongest(usym, self.symbol_count + self.opts.skip_amount);
            self.updateSCC(usym);
        }
        self.finalizeMean();
        self.finalizeOccurrences();
        self.finalizeLongest();
        self.finalizeChisq();
        self.finalizeEntropy();
        self.finalizeMonteCarlo();
        self.finalizeCompression();
        self.finalizeSCC();
        self.computeMarkov();
        if (self.opts.symbol_length != 8) { self.longest_byte_pos = (self.longest_position * self.opts.symbol_length) / 8; } else { self.longest_byte_pos = self.longest_position; }
    }

    // -- Output (writes directly to stdout) --
    pub fn printTerseOutput(self: *Self, terse_index: u32, filename: []const u8, parsed: ?fnparse.FilenameParsed) void {
    //pub fn printTerseOutput(self: *Self, terse_index: u32, parsed: ?fnparse.FilenameParsed) void {
        const r = self.results;
        if (self.opts.ent_exact) {
            const size = if (self.opts.symbol_length == 1) self.filebytes * 8 else self.filebytes;
            stdoutPrint("{d},{d},{d:.6},{d:.6},{d:.6},{d:.6},{d:.6}\n", .{ terse_index, size, r.entropy, r.chisq_distribution, r.mean, r.pi, r.scc });
        } else if (self.opts.parse_filename and self.opts.symbol_length == 1) {
            if (parsed) |p| {
                stdoutPrint("{d:>4},{d:>12},{s:>8},{s:>8},{d:>8.2},{d:>8.2},{d:>12.6},{d:>12.6},{x:>18},{d:>12.6},{d:>12.6},{d:>15.6},   {d:>16.6}, {d:>16.6}, {d:>16.6}, {d:>16.6}, {x:>19}, {d:>18}, {d:>18.6}, {d:>15}, {s}\n", .{ terse_index, self.symbol_count, p.device_id[0..p.device_id_len], p.process[0..p.process_len], p.voltage, p.temperature, r.entropy, r.min_entropy, r.min_entropy_symbol, r.chisq_percent, r.mean, r.pi, r.scc, r.p01, r.p10, self.markov_entropy, self.longest_longest_symbol, self.longest_longest, r.longest_pvalue, self.longest_byte_pos, filename });
            }
        } else {
            stdoutPrint("{d:>4},{d:>12},{d:>11.6}, {d:>12.6},{x:>18},{d:>12.6},{d:>12.6},{d:>15.6},       {d:>12.6}, {d:>16.6}, {d:>16.6}, {d:>16.6}, {x:>18}, {d:>18}, {d:>18.6}, {d:>15}, {s}\n", .{ terse_index, self.symbol_count, r.entropy, r.min_entropy, r.min_entropy_symbol, r.chisq_percent, r.mean, r.pi, r.scc, r.p01, r.p10, self.markov_entropy, self.longest_longest_symbol, self.longest_longest, r.longest_pvalue, self.longest_byte_pos, filename });
        }
    }

    //pub fn printVerboseOutput(self: *Self, filename: []const u8, parsed: ?fnparse.FilenameParsed) void {
    pub fn printVerboseOutput(self: *Self, parsed: ?fnparse.FilenameParsed) void {
        const r = self.results;
        const sl = self.opts.symbol_length;
        if (self.opts.print_occurrence) {
            for (0..@intCast(self.occurrence_size)) |i| { stdoutPrint("   Value {d:>4} , frequency={d} , fraction={d:.6}\n", .{ i, self.occurrence_count[i], @as(f64, @floatFromInt(self.occurrence_count[i])) / @as(f64, @floatFromInt(self.occurrence_total)) }); }
        }
        if (self.opts.print_longest) {
            for (0..@intCast(self.occurrence_size)) |i| { stdoutPrint("   Symbol {x} , Longest Run={d}\n", .{ i, self.longest_count[i] }); }
        }
        if (self.opts.ent_exact) {
            stdoutPrint("Entropy = {d:.6} bits per byte.\n\n", .{r.entropy * (8.0 / @as(f64, @floatFromInt(sl)))});
            stdoutPrint("Optimum compression would reduce the size\nof this {d} byte file by {d} percent\n\n", .{ self.filebytes, @as(i64, @intFromFloat(r.compression)) });
            stdoutPrint("Chi square distribution for {d} samples is {d:.6}, and randomly\nwould exceed this value less than {d:.6} percent of the times.\n\n", .{ r.chisq_count, r.chisq_distribution, r.chisq_percent });
            stdoutPrint("Arithmetic mean value of data bytes is {d:.6} (127.5 = random).\n", .{r.mean});
            stdoutPrint("Monte Carlo value for Pi is {d:.6} (error {d:.6} percent).\n", .{ r.pi, r.pi_err });
            stdoutPrint("Serial correlation coefficient is {d:.6} (totally uncorrelated = 0.0).\n", .{r.scc});
        } else {
            if (parsed) |p| { if (self.opts.parse_filename) { stdoutPrint("   Device ID   : {s}\n   Process     : {s}\n   Voltage     : {d:.2}V\n   Temperature : {d:.2}C\n", .{ p.device_id[0..p.device_id_len], p.process[0..p.process_len], p.voltage, p.temperature }); } }
            stdoutPrint("   Analysing {d} {d}-bit symbols\n", .{ self.symbol_count, sl });
            stdoutPrint("   Shannon IID Entropy = {d:.6} bits per symbol\n", .{r.entropy});
            stdoutPrint("   Optimal compression would compress by {d:.6} percent\n", .{r.compression});
            stdoutPrint("   Chi square: symbol count={d}, distribution={d:.2}, randomly exceeds {d:.2} percent of the time\n", .{ r.chisq_count, r.chisq_distribution, r.chisq_percent });
            stdoutPrint("   Mean = {d:.6}\n", .{r.mean});
            stdoutPrint("   Monte Carlo value for Pi is {d:.6} (error {d:.2} percent).\n", .{ r.pi, r.pi_err });
            stdoutPrint("   Serial Correlation = {d:.6}\n", .{r.scc});
            stdoutPrint("   Longest Run Symbol = {x}. Run Length = {d}\n", .{ self.longest_longest_symbol, self.longest_longest });
            if (sl == 1) stdoutPrint("   Probability of longest run being <= {d} = {d:.6}\n", .{ self.longest_longest, r.longest_pvalue });
            stdoutPrint("   Position of Longest Run = {d} (0x{x}). Byte position {d} (0x{x})\n", .{ self.longest_position, self.longest_position, self.longest_byte_pos, self.longest_byte_pos });
            stdoutPrint("   A 2 state Markov generator with transition probabilities P01={d:.6}, P10={d:.6} would generate data with entropy {d:.6} per bit with 8 bit symbols with the same mean and serial correlation\n", .{ r.p01, r.p10, self.markov_entropy });
        }
    }
};
