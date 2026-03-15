// Converted from C++ to Rust

use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;

pub fn hash_combine<T: Hash>(seed: &mut u64, value: T) {
    let mut hasher = DefaultHasher::new();
    value.hash(&mut hasher);
    let hash = hasher.finish();
    *seed ^= hash.wrapping_add(0x9e3779b9)
        .wrapping_add(*seed << 6)
        .wrapping_add(*seed >> 2);
}

pub fn hash_range<T: Hash>(items: &[T]) -> u64 {
    let mut seed = 0u64;
    for item in items {
        hash_combine(&mut seed, item);
    }
    seed
}
