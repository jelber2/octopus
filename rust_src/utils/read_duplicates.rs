// Converted from C++ to Rust

use crate::basics::aligned_read::AlignedRead;

pub fn mark_duplicates(reads: &mut Vec<AlignedRead>) {
    // Placeholder: real implementation would use read coordinates and sequences
    let _ = reads;
}

pub fn remove_duplicates(reads: Vec<AlignedRead>) -> Vec<AlignedRead> {
    let mut seen = std::collections::HashSet::new();
    reads.into_iter().filter(|r| {
        let key = (r.name().to_string(), r.region().begin(), r.region().end());
        seen.insert(key)
    }).collect()
}
