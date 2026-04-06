const std = @import("std");

const alloc = std.heap.smp_allocator.alloc;


const Vec3 = [3]f64;

const Mesh = struct
{
    vertices: []Vec3,
    indices:  []u32,    

    pub fn init (vertices_len:usize, indices_len:usize) !@This() {
        return .{
            .vertices = try alloc(Vec3, vertices_len),
            .indices  = try alloc(u32, indices_len),
        };
    }
    pub fn Mesh(vertices_len:usize, indices_len:usize) !@This() {
        return .{
            .vertices = try alloc(Vec3, vertices_len),
            .indices  = try alloc(u32, indices_len),
        };
    }
};

fn sum (x_values:[]f64) f64
{
    var r: f64 = 0;
    for (x_values) |x|
    {
        r += x;
    }

    return r;
}


fn main () !void
{
    const len = 100;
    // const size = len * @sizeOf(f64);
    var buf = try std.heap.smp_allocator.alloc(f64, len);
    @memset(buf, 0);

    const mesh: Mesh = try .init(1000, 3000);
    // const mesh: Mesh = try .Mesh(1000, 3000);
    
    const vertices_len = mesh.vertices.len;
    const indices_count = mesh.indices.len / 3;

    for (0..vertices_len) |i| {
        mesh.vertices[i] = .{1,1,1};
    }

    for (0..indices_count) |i| {
        mesh.indices[3*i + 0] = i + 1;
        mesh.indices[3*i + 1] = i + 2;
        mesh.indices[3*i + 2] = i + 3;
    }     

    for (0..len) |i| {
        const x: f64 = @floatFromInt(i);
        buf[i] = (0.2*x*x + 2*x);
    }

    const r = sum(buf);
}
