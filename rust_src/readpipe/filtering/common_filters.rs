// Converted from C++ to Rust

use crate::basics::aligned_read::AlignedRead;
use super::read_filter::ReadFilter;

pub struct MinMappingQualityFilter {
    min_quality: u8,
}

impl MinMappingQualityFilter {
    pub fn new(min_quality: u8) -> Self { MinMappingQualityFilter { min_quality } }
}

impl ReadFilter for MinMappingQualityFilter {
    fn name(&self) -> &str { "MinMappingQuality" }
    fn passes(&self, read: &AlignedRead) -> bool {
        read.mapping_quality() >= self.min_quality
    }
}

pub struct UnmappedReadFilter;

impl ReadFilter for UnmappedReadFilter {
    fn name(&self) -> &str { "UnmappedRead" }
    fn passes(&self, read: &AlignedRead) -> bool {
        !read.is_unmapped()
    }
}

pub struct DuplicateReadFilter;

impl ReadFilter for DuplicateReadFilter {
    fn name(&self) -> &str { "DuplicateRead" }
    fn passes(&self, read: &AlignedRead) -> bool {
        !read.is_duplicate()
    }
}

pub struct SecondaryAlignmentFilter;

impl ReadFilter for SecondaryAlignmentFilter {
    fn name(&self) -> &str { "SecondaryAlignment" }
    fn passes(&self, read: &AlignedRead) -> bool {
        !read.is_secondary_alignment()
    }
}

pub struct SupplementaryAlignmentFilter;

impl ReadFilter for SupplementaryAlignmentFilter {
    fn name(&self) -> &str { "SupplementaryAlignment" }
    fn passes(&self, read: &AlignedRead) -> bool {
        !read.is_supplementary_alignment()
    }
}

pub struct MinReadLengthFilter {
    min_length: usize,
}

impl MinReadLengthFilter {
    pub fn new(min_length: usize) -> Self { MinReadLengthFilter { min_length } }
}

impl ReadFilter for MinReadLengthFilter {
    fn name(&self) -> &str { "MinReadLength" }
    fn passes(&self, read: &AlignedRead) -> bool {
        read.sequence().len() >= self.min_length
    }
}
