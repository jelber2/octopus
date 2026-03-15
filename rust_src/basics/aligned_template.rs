// Converted from C++ to Rust

use super::aligned_read::AlignedRead;

#[derive(Debug, Clone)]
pub struct AlignedTemplate {
    reads: Vec<AlignedRead>,
}

impl AlignedTemplate {
    pub fn new(reads: Vec<AlignedRead>) -> Self {
        AlignedTemplate { reads }
    }

    pub fn reads(&self) -> &[AlignedRead] { &self.reads }
    pub fn len(&self) -> usize { self.reads.len() }
    pub fn is_empty(&self) -> bool { self.reads.is_empty() }
}
