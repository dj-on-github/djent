//! djent - A reimplementation of Fourmilab's ent with several improvements.
//!
//! Copyright (C) 2017  David Johnston
//! Translated to Zig 0.15.x, 2026
//!
//! Contact. David Johnston dj@deadhat.com

const std = @import("std");
const Analyzer = @import("analyzer.zig").Analyzer;
const Options = @import("analyzer.zig").Options;
const fnparse = @import("filename_parse.zig");

// Simple stdout helper for Zig 0.15.
fn stdoutPrint(comptime fmt: []const u8, args: anytype) void {
    var buf: [8192]u8 = undefined;
    const msg = std.fmt.bufPrint(&buf, fmt, args) catch return;
    const f = std.fs.File.stdout();
    f.writeAll(msg) catch {};
}

fn displayUsage() void {
    std.debug.print(
        \\Usage: djent [-brRpcCuhds] [-l <n>] [-i <input file list filename>] [filename] [filename2] ...
        \\
        \\Compute statistics of random data.
        \\  Author: David Johnston, dj@deadhat.com
        \\
        \\  -i <filename>  --inputfilelist=<filename> Read list of filenames from <filename>
        \\  -p             --parse_filename           Extract CID, Process, Voltage and Temperature from filename.
        \\  -l <n>         --symbol_length=<n>        Treat incoming data symbols as bitlength n. Default is 8.
        \\  -b             --binary                   Treat incoming data as binary. Default bit length will be -l 1
        \\  -r             --byte_reverse             Reverse the bit order in incoming bytes
        \\  -R             --word_reverse             Reverse the byte order in incoming 4 byte words
        \\  -c             --occurrence               Print symbol occurrence counts
        \\  -C             --longest                  Print symbol longest run counts
        \\  -w             --scc_wrap                 Treat data as cyclical in SCC
        \\  -n <n>         --lagn=<n>                 Lag gap in SCC. Default=1
        \\  -S <n>         --skip=<n>                 Skip over <n> initial symbols
        \\  -L <n>         --substring=<n>            Analyse no more than <n> symbols
        \\  -f             --fold                     Fold uppercase letters to lower case
        \\  -t             --terse                    Terse output
        \\  -e             --ent_exact                Exactly match output format of ent
        \\  -s             --suppress_header          Suppress the header in terse output
        \\  -h or -u       --help                     Print this text
        \\
    , .{});
}

fn printTerseHeader(opts: Options) void {
    if (opts.ent_exact) {
        if (opts.symbol_length == 1) { stdoutPrint("0,File-bits,Entropy,Chi-square,Mean,Monte-Carlo-Pi,Serial-Correlation\n", .{}); } else { stdoutPrint("0,File-bytes,Entropy,Chi-square,Mean,Monte-Carlo-Pi,Serial-Correlation\n", .{}); }
    } else if (opts.parse_filename) {
        stdoutPrint("   0,     symbols,     CID, Process, Voltage,    Temp,     Entropy,  MinEntropy, MinEntropy-Symbol,  Chi-square,        Mean, Monte-Carlo-Pi, Serial-Correlation,              P01,              P10,  mkv_min_entropy,  Longest-Run-Symbol, Longest-Run-Length, Longest-Run-PValue, Longest-Run-Pos, Filename\n", .{});
    } else {
        stdoutPrint("   0,     symbols,    Entropy,  Min_entropy, MinEntropy-Symbol,  Chi-square,        Mean, Monte-Carlo-Pi, Serial-Correlation,              P01,              P10,  mkv_min_entropy, Longest-Run-Symbol, Longest-Run-Length, Longest-Run-PValue, Longest-Run-Pos, Filename\n", .{});
    }
}

/// Read a file list from a text file, skipping lines starting with #.
fn readFileList(allocator: std.mem.Allocator, filename: []const u8) !std.ArrayListUnmanaged([]const u8) {
    var list: std.ArrayListUnmanaged([]const u8) = .{};
    const data = try std.fs.cwd().readFileAlloc(allocator, filename, 10 * 1024 * 1024);
    defer allocator.free(data);

    var lines = std.mem.splitScalar(u8, data, '\n');
    while (lines.next()) |line| {
        const trimmed = std.mem.trim(u8, line, &[_]u8{ ' ', '\t', '\r' });
        if (trimmed.len == 0 or trimmed[0] == '#') continue;
        try list.append(allocator, try allocator.dupe(u8, trimmed));
    }
    return list;
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    // Parse command-line arguments
    var opts = Options{};
    var filenames: std.ArrayListUnmanaged([]const u8) = .{};
    defer {
        for (filenames.items) |f| allocator.free(f);
        filenames.deinit(allocator);
    }
    var got_symbol_length = false;

    var args = std.process.args();
    _ = args.next(); // skip program name

    while (args.next()) |arg| {
        if (arg.len > 1 and arg[0] == '-' and arg[1] != '-') {
            // Short options
            var ci: usize = 1;
            while (ci < arg.len) {
                const c = arg[ci];
                switch (c) {
                    'b' => { if (!got_symbol_length) opts.symbol_length = 1; opts.hex_mode = false; },
                    'l' => { if (args.next()) |v| { opts.symbol_length = std.fmt.parseInt(u32, v, 10) catch 8; got_symbol_length = true; } },
                    'i' => { if (args.next()) |v| { opts.input_list_filename = try allocator.dupe(u8, v); opts.using_input_list_file = true; } },
                    'c' => opts.print_occurrence = true,
                    'C' => opts.print_longest = true,
                    'p' => opts.parse_filename = true,
                    'r' => opts.byte_reverse = true,
                    'R' => opts.word_reverse = true,
                    'w' => opts.scc_wrap = true,
                    'n' => { if (args.next()) |v| { opts.lag_n = std.fmt.parseInt(u32, v, 10) catch 1; } },
                    'f' => opts.fold = true,
                    't' => opts.terse = true,
                    'e' => opts.ent_exact = true,
                    's' => opts.suppress_header = true,
                    'S' => { if (args.next()) |v| { opts.got_skip = true; opts.skip_amount = std.fmt.parseInt(u64, v, 10) catch 0; } },
                    'L' => { if (args.next()) |v| { opts.got_substring = true; opts.substring = std.fmt.parseInt(u64, v, 10) catch 0; } },
                    'h', 'u', '?' => { displayUsage(); std.process.exit(0); },
                    else => {},
                }
                ci += 1;
            }
        } else if (arg.len > 2 and std.mem.startsWith(u8, arg, "--")) {
            const option = arg[2..];
            if (std.mem.eql(u8, option, "help")) {
                displayUsage(); std.process.exit(0); 
            } else if (std.mem.eql(u8, option, "binary")) {
                if (!got_symbol_length) opts.symbol_length = 1;
                opts.hex_mode = false; 
            } else if (std.mem.startsWith(u8, option, "symbol_length=")) {
                opts.symbol_length = std.fmt.parseInt(u32, option["symbol_length=".len..], 10) catch 8;
                got_symbol_length = true;
            } else if (std.mem.startsWith(u8, option, "inputfilelist=")) {
                opts.input_list_filename = try allocator.dupe(u8, option["inputfilelist=".len..]);
                opts.using_input_list_file = true; 
            } else if (std.mem.eql(u8, option, "occurrence")) { 
                opts.print_occurrence = true; 
            } else if (std.mem.eql(u8, option, "longest")) { 
                opts.print_longest = true; 
            } else if (std.mem.eql(u8, option, "parse_filename")) { 
                opts.parse_filename = true; 
            } else if (std.mem.eql(u8, option, "byte_reverse")) { 
                opts.byte_reverse = true; 
            } else if (std.mem.eql(u8, option, "word_reverse")) { 
                opts.word_reverse = true; 
            } else if (std.mem.eql(u8, option, "scc_wrap")) { 
                opts.scc_wrap = true; 
            } else if (std.mem.startsWith(u8, option, "lagn=")) { 
                opts.lag_n = std.fmt.parseInt(u32, option["lagn=".len..], 10) catch 1; 
            } else if (std.mem.eql(u8, option, "fold")) { 
                opts.fold = true; 
            } else if (std.mem.eql(u8, option, "terse")) { 
                opts.terse = true; 
            } else if (std.mem.eql(u8, option, "ent_exact")) { 
                opts.ent_exact = true; 
            } else if (std.mem.eql(u8, option, "suppress_header")) { 
                opts.suppress_header = true; 
            } else if (std.mem.startsWith(u8, option, "skip=")) { 
                opts.got_skip = true; 
                opts.skip_amount = std.fmt.parseInt(u64, option["skip=".len..], 10) catch 0; 
            } else if (std.mem.startsWith(u8, option, "substring=")) { 
                opts.got_substring = true; 
                opts.substring = std.fmt.parseInt(u64, option["substring=".len..], 10) catch 0; 
            }
        } else {
            try filenames.append(allocator, try allocator.dupe(u8, arg));
        }
    }

    opts.use_stdin = (filenames.items.len == 0 and !opts.using_input_list_file);

    // Validation
    if (opts.fold and opts.symbol_length != 8) { std.debug.print("Error: Fold must be used with 8 bit word size\n", .{}); std.process.exit(1); }
    if (opts.symbol_length < 1) { std.debug.print("Error: Symbol length must not be 0.\n", .{}); std.process.exit(1); }
    if (opts.parse_filename and opts.use_stdin) { std.debug.print("Error: Can't parse filename when using stdin\n", .{}); std.process.exit(1); }
    if (opts.got_skip and opts.skip_amount < 1) { std.debug.print("Error: skip amount must be > 0\n", .{}); std.process.exit(1); }
    if (opts.got_substring and opts.substring < 1) { std.debug.print("Error: substring length must be > 0\n", .{}); std.process.exit(1); }

    // Build file list
    var file_list: std.ArrayListUnmanaged([]const u8) = .{};
    var own_file_list = false;
    defer {
        if (own_file_list) {
            for (file_list.items) |item| allocator.free(item);
            file_list.deinit(allocator);
        }
    }

    if (opts.using_input_list_file) {
        file_list = readFileList(allocator, opts.input_list_filename) catch {
            std.debug.print("Error: Cannot open {s} for reading\n", .{opts.input_list_filename});
            std.process.exit(1);
        };
        own_file_list = true;
        if (file_list.items.len == 0) { std.debug.print("Error: No filenames found in {s}\n", .{opts.input_list_filename}); std.process.exit(1); }
    } else if (opts.use_stdin) {
        try file_list.append(allocator, try allocator.dupe(u8, "<stdin>"));
        own_file_list = true;
    } else {
        // Use filenames directly (already owned by the filenames list)
        file_list = filenames;
        own_file_list = false;
    }

    var terse_index: u32 = 0;

    for (file_list.items) |filename| {
        terse_index += 1;

        var az = Analyzer.create(allocator, opts);
        defer az.deinit();

        if (std.mem.eql(u8, filename, "<stdin>")) {
            if (opts.terse and terse_index == 1 and !opts.suppress_header) printTerseHeader(opts);
            const stdin_file = std.fs.File.stdin();
            try az.process(stdin_file);
            if (opts.terse) { az.printTerseOutput(terse_index, filename, null); } else { az.printVerboseOutput(null); }
        } else {
            if (opts.terse and terse_index == 1 and !opts.suppress_header) printTerseHeader(opts);
            if (!opts.terse and !opts.ent_exact) {
                const mode = if (opts.hex_mode) "hex text" else "binary";
                stdoutPrint(" opening {s} as {s}\n Symbol Size(bits) = {d}\n", .{ filename, mode, opts.symbol_length });
            }

            const file = std.fs.cwd().openFile(filename, .{}) catch {
                std.debug.print("Error: Unable to open file {s}\n", .{filename});
                std.process.exit(1);
            };
            defer file.close();

            const parsed_fn: ?fnparse.FilenameParsed = if (opts.parse_filename) fnparse.parseFilename(filename) else null;

            try az.process(file);
            if (opts.terse) { az.printTerseOutput(terse_index, filename, parsed_fn); } else { az.printVerboseOutput( parsed_fn); }
        }
    }
}
