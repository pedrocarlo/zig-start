//! Implementation as seen in Python std lib with some more comments, but in zig.
//! However, not all convinience methods are implemented
//! https://github.com/python/cpython/blob/3.13/Lib/heapq.py

const std = @import("std");
const Allocator = std.mem.Allocator;
const assert = std.debug.assert;
const testing = std.testing;
const random = std.rand;
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
        compareFn: @TypeOf(compareFn),

        /// Deinitialize with `deinit`
        pub fn init(allocator: Allocator) Self {
            return Self{ .array = std.ArrayList(T).init(allocator), .allocator = allocator, .compareFn = compareFn };
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
        fn heapify(self: Self) void {
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
        pub fn pushSlice(self: Self, items: []const T) Allocator.Error!void {
            try self.array.appendSlice(items);
            // After appending adjust the heap
            self.heapify();
        }

        /// Push element into heap
        pub fn push(self: Self, val: T) Allocator.Error!void {
            try self.array.append(val);
            self.siftdown(0, self.size());
        }

        /// Pop element from heap.
        ///
        /// Returns optional T as you can pop the heap when its empty
        pub fn pop(self: Self) ?T {
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
        pub fn peek(self: Self) ?T {
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

        fn siftdown(self: Self, start_pos: usize, pos: usize) void {
            const new_item = self.get(pos);

            // While we are lower in the tree than the node in start_pos
            while (pos > start_pos) {
                // Same as Floor((pos-1)/2) to get the parent
                const parent_pos: usize = (pos - 1) >> 1;
                const parent = self.get(parent_pos);

                // If new item is smaller than parent, swap parent pos
                if (compareFn(new_item, parent) == .lt) {
                    self.set(pos, parent);
                    pos = parent_pos;
                    continue;
                }
                break;
            }
            // Here we got top most position we needed to get to
            self.set(pos, new_item);
        }

        fn siftup(self: Self, pos: usize) void {
            const end_pos = self.size();
            const start_pos = pos;
            const new_item = self.get(pos);
            // Bubble up the smaller child until hitting a leaf.
            var child_pos = 2 * pos + 1; // leftmost child position

            // while we are at higher in the tree than the last leaf
            while (child_pos > end_pos) {
                const right_pos = child_pos + 1;

                // If right_pos has a parent array[right_pos] is smaller than value array[child_pos]
                // make right_pos to be the smallest value index to follow
                if (right_pos < end_pos and !self.get(child_pos) < self.get(right_pos)) {
                    child_pos = right_pos;
                }
                // Move the smaller child up
                self.set(pos, self.get(child_pos));
                pos = child_pos;
                // Get left child
                child_pos = 2 * pos + 1;
            }
            // The leaf at pos is empty now. Put newitem there, and bubble it up
            // to its final resting place (by sifting its parents down).
            self.set(pos, new_item);

            self.siftdown(start_pos, pos);
        }

        /// Get value at index
        fn get(self: Self, index: usize) T {
            assert(self.array.items.len > index);
            return self.array.items[index] orelse unreachable;
        }

        /// Set value at index
        fn set(self: Self, index: usize, val: T) void {
            assert(self.array.items.len > index);
            self.array.items[index] = val;
        }
    };
}

// Tests gotten from Python std lib as well

fn compare(a: u32, b: u32) Order {
    return std.math.order(a, b);
}

const MinHeap = BinaryHeapArrayList(u32, compare);

/// Testing helper to check if the heap continues to be a heap
fn checkInvariant(comptime T: type, comptime compareFn: fn (a: T, b: T) Order, Heap: BinaryHeapArrayList(T, compareFn)) !void {
    for (Heap.array.items, 0..) |item, index| {
        if (index > 0) {
            const parent_pos = (index - 1) >> 1;
            const comparison = Heap.compareFn(Heap.get(parent_pos), item);
            try testing.expect(comparison.compare(.lte));
        }
    }
}

test "push pop" {
    const heap = MinHeap.init(testing.allocator);
    defer heap.deinit();
    try checkInvariant(u32, compare, heap);
}
