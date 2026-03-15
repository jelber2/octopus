// Converted from C++ to Rust

use std::collections::HashMap;

pub struct KmerMapper {
    k: usize,
    map: HashMap<Vec<u8>, Vec<u32>>,
}

impl KmerMapper {
    pub fn new(k: usize) -> Self {
        KmerMapper { k, map: HashMap::new() }
    }

    pub fn add(&mut self, sequence: &[u8], position: u32) {
        for i in 0..sequence.len().saturating_sub(self.k - 1) {
            let kmer = sequence[i..i + self.k].to_vec();
            self.map.entry(kmer).or_default().push(position + i as u32);
        }
    }

    pub fn find(&self, kmer: &[u8]) -> Option<&[u32]> {
        self.map.get(kmer).map(|v| v.as_slice())
    }

    pub fn k(&self) -> usize { self.k }
    pub fn size(&self) -> usize { self.map.len() }
}
