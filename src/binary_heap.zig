//! Implementation as seen in Python std lib with some more comments, but in zig.
//! However, not all convinience methods are implemented
//! https://github.com/python/cpython/blob/3.13/Lib/heapq.py

const std = @import("std");
const Allocator = std.mem.Allocator;
const assert = std.debug.assert;
const testing = std.testing;
const Order = std.math.Order;

/// Binary Min Heap Implementation with an underlying
/// Array List data structure.
///
/// `compareFn` must return Order type from stdlib.
pub fn BinaryHeapArrayList(comptime T: type, comptime compareFn: fn (a: T, b: T) Order) type {
    return struct {
        const Self = @This();

        /// Underlying array list
        array: std.ArrayList(T),
        allocator: Allocator,

        /// Deinitialize with `deinit`
        pub fn init(allocator: Allocator) Self {
            return Self{ .array = std.ArrayList(T).init(allocator), .allocator = allocator };
        }

        /// Release all allocated memory
        pub fn deinit(self: Self) void {
            // Got this from std
            // Probably used to safeguard against storing some null type
            if (@sizeOf(T) > 0) {
                // Array List manages the memory
                self.array.deinit();
            }
        }

        /// From python docs:
        ///
        /// Heapify current heap, in-place, in O(len(heap)) time.
        fn heapify(self: *Self) void {
            // From python docs:
            // Transform bottom-up. The largest index there's any point to looking at
            // is the largest with a child index in-range, so must have 2*i + 1 < n,
            // or i < (n-1)/2. If n is even = 2*j, this is (2*j-1)/2 = j-1/2 so
            // j-1 is the largest, which is n//2 - 1. If n is odd = 2*j+1, this is
            // (2*j+1-1)/2 = j so j-1 is the largest, and that's again n//2-1.
            const n = self.size() >> 1; // same as n // 2
            while (n > 0) {
                n -= 1;
                self.siftup(n);
            }
        }

        /// Push many elements. Heap is subsequently heapified to maintain properties
        pub fn pushSlice(self: *Self, items: []const T) Allocator.Error!void {
            try self.array.appendSlice(items);
            // After appending adjust the heap
            self.heapify();
        }

        /// Push element into heap
        pub fn push(self: *Self, val: T) Allocator.Error!void {
            try self.array.append(val);
            self.siftdown(0, self.size() - 1);
        }

        /// Pop element from heap.
        ///
        /// Returns T
        pub fn pop(self: *Self) T {
            // This will be the value to return if the heap is only 1 item
            const last = self.array.pop();
            if (self.size() > 0) {
                // Get the head of heap
                const return_item = self.get(0);
                // Set last value as head and sift up the heap
                // Last will be placed in the correct spot as with
                // siftup implementaion it calls siftdown on what was displaced
                self.set(0, last);
                self.siftup(0);
                return return_item;
            }
            return last;
        }

        /// Pop element from heap.
        ///
        /// Returns optional T as you can pop the heap when its empty
        pub fn popOrNull(self: *Self) ?T {
            // This will be the value to return if the heap is only 1 item
            const last = self.array.popOrNull();
            if (self.size() > 0) {
                // Get the head of heap
                const return_item = self.get(0);
                // Set last value as head and sift up the heap
                // Last will be placed in the correct spot as with
                // siftup implementaion it calls siftdown on what was displaced
                self.set(0, last);
                self.siftup(0);
                return return_item;
            }
            return last;
        }

        /// Peeks at the head of the heap
        pub fn peek(self: *Self) ?T {
            if (self.size() == 0) {
                return null;
            }
            return self.get(0);
        }

        /// Size of heap
        pub fn size(self: Self) usize {
            return self.array.items.len;
        }

        /// Returns a Clone of the underlying `std.ArrayList`
        pub fn toArrayList(self: Self) Allocator.Error!std.ArrayList(T) {
            return try self.array.clone();
        }

        fn siftdown(self: *Self, start_pos: usize, pos: usize) void {
            const new_item = self.get(pos);

            var func_pos = pos;

            // While we are lower in the tree than the node in start_pos
            while (func_pos > start_pos) {
                // Same as Floor((pos-1)/2) to get the parent
                const parent_pos: usize = (func_pos - 1) >> 1;
                const parent = self.get(parent_pos);

                // If new item is smaller than parent, swap parent pos
                if (compareFn(new_item, parent) == .lt) {
                    self.set(func_pos, parent);
                    func_pos = parent_pos;
                    continue;
                }
                break;
            }
            // Here we got top most position we needed to get to
            self.set(func_pos, new_item);
        }

        /// Swaps item down starting from pos
        fn siftup(self: *Self, pos: usize) void {
            const end_pos = self.size();
            const start_pos = pos;
            var func_pos = pos;
            const new_item = self.get(func_pos);
            // Bubble up the smaller child until hitting a leaf.
            var child_pos = (2 * func_pos) + 1; // leftmost child position

            // while we are lower in the tree than the end_pos
            while (child_pos < end_pos) {
                const right_pos = child_pos + 1;

                // If right_pos has a parent array[right_pos] is smaller than value array[child_pos]
                // make right_pos to be the smallest value index to follow
                if (right_pos < end_pos and !compareFn(self.get(child_pos), self.get(right_pos)).compare(.lt)) {
                    child_pos = right_pos;
                }
                // Move the smaller child up
                self.set(func_pos, self.get(child_pos));
                func_pos = child_pos;
                // Get left child
                child_pos = (2 * func_pos) + 1;
            }
            // The leaf at pos is empty now. Put newitem there, and bubble it up
            // to its final resting place (by sifting its parents down).
            self.set(func_pos, new_item);

            self.siftdown(start_pos, func_pos);
        }

        /// Get value at index
        fn get(self: Self, index: usize) T {
            // std.debug.print("\nlen: {d}, index: {d}\n", .{ self.array.items.len, index });
            assert(self.array.items.len > index);
            return self.array.items[index];
        }

        /// Set value at index
        fn set(self: *Self, index: usize, val: T) void {
            assert(self.array.items.len > index);
            self.array.items[index] = val;
        }
    };
}

pub const std_options = .{
    // Set the log level to info
    .log_level = .info,
};

// Tests gotten from Python std lib as well

fn compare(a: u32, b: u32) Order {
    return std.math.order(a, b);
}

const MinHeap = BinaryHeapArrayList(u32, compare);

/// Testing helper to check if the heap continues to be a heap
fn checkInvariant(comptime T: type, comptime compareFn: fn (a: T, b: T) Order, Heap: []T) !void {
    for (Heap, 0..) |item, index| {
        if (index > 0) {
            const parent_pos = (index - 1) >> 1;
            const comparison = compareFn(Heap[parent_pos], item);
            try testing.expect(comparison.compare(.lte));
        }
    }
}

fn lessThan(context: @TypeOf(.{}), a: u32, b: u32) bool {
    _ = context;
    return a < b;
}

test "simple push pop" {
    const a = testing.allocator;
    var heap = MinHeap.init(a);
    defer heap.deinit();

    try heap.push(1);
    try checkInvariant(u32, compare, heap.array.items);

    try heap.push(2);
    try checkInvariant(u32, compare, heap.array.items);

    try heap.push(3);
    try checkInvariant(u32, compare, heap.array.items);

    try testing.expect(heap.pop() == 1);
    // std.debug.print("{d}\n", .{heap.array.items[0..]});
    try checkInvariant(u32, compare, heap.array.items);

    try testing.expect(heap.pop() == 2);
    try checkInvariant(u32, compare, heap.array.items);

    try testing.expect(heap.pop() == 3);
    try checkInvariant(u32, compare, heap.array.items);
}

test "push pop" {
    const a = testing.allocator;
    var heap = MinHeap.init(a);
    defer heap.deinit();
    try checkInvariant(u32, compare, heap.array.items);
    var prng = std.rand.DefaultPrng.init(blk: {
        var seed: u64 = undefined;
        try std.posix.getrandom(std.mem.asBytes(&seed));
        break :blk seed;
    });
    var data = std.mem.zeroes([256]u32);
    const rand = prng.random();
    for (0..256) |i| {
        const item = rand.int(u32);

        data[i] = item;
        try heap.push(item);
        try checkInvariant(u32, compare, heap.array.items);
    }
    var results = std.ArrayList(u32).init(a);
    defer results.deinit();

    while (heap.size() > 0) {
        // std.debug.print("heap size: {d}\n", .{heap.size()});
        const item = heap.pop();
        try checkInvariant(u32, compare, heap.array.items);
        try results.append(item);
    }

    std.sort.heap(u32, &data, .{}, lessThan);
    // std.debug.print("data :{d}\n", .{data[0..]});
    // std.debug.print("result: {d}\n", .{results.items[0..]});

    try testing.expectEqualDeep(results.items, &data);

    try checkInvariant(u32, compare, results.items);
}
