const std = @import("std");

pub fn build(b: *std.Build) void {
    const exe = b.addExecutable(.{
        .name = "djent",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = b.standardTargetOptions(.{}),
            .optimize = b.standardOptimizeOption(.{}),
        }),
    });

    b.installArtifact(exe);

    const run_cmd = b.addRunArtifact(exe);
    
    const run_step = b.step("run", "Run djent");

    run_step.dependOn(&run_cmd.step);
}
