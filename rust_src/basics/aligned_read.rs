// Converted from C++ to Rust

use super::genomic_region::GenomicRegion;
use super::cigar_string::CigarString;

pub type MappingQuality = u8;
pub type BaseQuality = u8;
pub type Sequence = Vec<u8>;

#[derive(Debug, Clone)]
pub struct AlignedRead {
    name: String,
    region: GenomicRegion,
    sequence: Sequence,
    qualities: Vec<BaseQuality>,
    cigar: CigarString,
    mapping_quality: MappingQuality,
    flags: u16,
}

impl AlignedRead {
    pub fn new(
        name: String,
        region: GenomicRegion,
        sequence: Sequence,
        qualities: Vec<BaseQuality>,
        cigar: CigarString,
        mapping_quality: MappingQuality,
        flags: u16,
    ) -> Self {
        AlignedRead { name, region, sequence, qualities, cigar, mapping_quality, flags }
    }

    pub fn name(&self) -> &str { &self.name }
    pub fn region(&self) -> &GenomicRegion { &self.region }
    pub fn mapped_region(&self) -> &GenomicRegion { &self.region }
    pub fn sequence(&self) -> &Sequence { &self.sequence }
    pub fn qualities(&self) -> &[BaseQuality] { &self.qualities }
    pub fn qualities_mut(&mut self) -> &mut Vec<BaseQuality> { &mut self.qualities }
    pub fn cigar(&self) -> &CigarString { &self.cigar }
    pub fn mapping_quality(&self) -> MappingQuality { self.mapping_quality }
    pub fn flags(&self) -> u16 { self.flags }

    pub fn is_paired(&self) -> bool { self.flags & 0x1 != 0 }
    pub fn is_proper_pair(&self) -> bool { self.flags & 0x2 != 0 }
    pub fn is_unmapped(&self) -> bool { self.flags & 0x4 != 0 }
    pub fn mate_is_unmapped(&self) -> bool { self.flags & 0x8 != 0 }
    pub fn is_reverse_strand(&self) -> bool { self.flags & 0x10 != 0 }
    pub fn mate_is_reverse_strand(&self) -> bool { self.flags & 0x20 != 0 }
    pub fn is_first_in_pair(&self) -> bool { self.flags & 0x40 != 0 }
    pub fn is_second_in_pair(&self) -> bool { self.flags & 0x80 != 0 }
    pub fn is_secondary_alignment(&self) -> bool { self.flags & 0x100 != 0 }
    pub fn is_qc_fail(&self) -> bool { self.flags & 0x200 != 0 }
    pub fn is_duplicate(&self) -> bool { self.flags & 0x400 != 0 }
    pub fn is_supplementary_alignment(&self) -> bool { self.flags & 0x800 != 0 }
}
