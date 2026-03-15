// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Converted from C++ to Rust

use std::fmt;
use std::hash::{Hash, Hasher};
use crate::basics::genomic_region::GenomicRegion;
use crate::basics::contig_region::ContigRegion;

pub type NucleotideSequence = Vec<u8>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContigAllele {
    pub region: ContigRegion,
    pub sequence: NucleotideSequence,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Allele {
    pub region: GenomicRegion,
    pub sequence: NucleotideSequence,
}

impl ContigAllele {
    pub fn new(region: ContigRegion, sequence: NucleotideSequence) -> Self {
        ContigAllele { region, sequence }
    }

    pub fn from_position(begin: u32, sequence: NucleotideSequence) -> Result<Self, crate::basics::contig_region::BadRegion> {
        let end = begin + sequence.len() as u32;
        let region = ContigRegion::new(begin, end)?;
        Ok(ContigAllele { region, sequence })
    }

    pub fn mapped_region(&self) -> &ContigRegion { &self.region }
    pub fn sequence(&self) -> &NucleotideSequence { &self.sequence }
    pub fn sequence_size(&self) -> usize { self.sequence.len() }
    pub fn is_sequence_empty(&self) -> bool { self.sequence.is_empty() }
}

impl Allele {
    pub fn new(region: GenomicRegion, sequence: NucleotideSequence) -> Self {
        Allele { region, sequence }
    }

    pub fn from_parts(contig: impl Into<String>, begin: u32, sequence: NucleotideSequence) -> Result<Self, String> {
        let end = begin + sequence.len() as u32;
        let region = GenomicRegion::new(contig, begin, end).map_err(|e| e.to_string())?;
        Ok(Allele { region, sequence })
    }

    pub fn mapped_region(&self) -> &GenomicRegion { &self.region }
    pub fn sequence(&self) -> &NucleotideSequence { &self.sequence }
    pub fn sequence_size(&self) -> usize { self.sequence.len() }
    pub fn is_sequence_empty(&self) -> bool { self.sequence.is_empty() }
}

impl PartialOrd for ContigAllele {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ContigAllele {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.region.cmp(&other.region)
            .then_with(|| self.sequence.cmp(&other.sequence))
    }
}

impl Hash for ContigAllele {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.region.hash(state);
        self.sequence.hash(state);
    }
}

impl PartialOrd for Allele {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.region.contig_name() != other.region.contig_name() {
            return None;
        }
        Some(self.region.contig_region().cmp(other.region.contig_region())
            .then_with(|| self.sequence.cmp(&other.sequence)))
    }
}

impl Hash for Allele {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.region.hash(state);
        self.sequence.hash(state);
    }
}

impl fmt::Display for ContigAllele {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.region, String::from_utf8_lossy(&self.sequence))
    }
}

impl fmt::Display for Allele {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.region, String::from_utf8_lossy(&self.sequence))
    }
}

pub fn is_insertion_allele(allele: &Allele) -> bool {
    allele.sequence.len() > crate::basics::contig_region::size(allele.region.contig_region()) as usize
}

pub fn is_deletion_allele(allele: &Allele) -> bool {
    allele.sequence.len() < crate::basics::contig_region::size(allele.region.contig_region()) as usize
}

pub fn is_indel_allele(allele: &Allele) -> bool {
    is_insertion_allele(allele) || is_deletion_allele(allele)
}

pub fn is_simple_insertion(allele: &Allele) -> bool {
    is_insertion_allele(allele) && crate::basics::contig_region::is_empty(allele.region.contig_region())
}

pub fn is_simple_deletion(allele: &Allele) -> bool {
    is_deletion_allele(allele) && allele.sequence.is_empty()
}

pub fn demote(allele: Allele) -> ContigAllele {
    ContigAllele { region: *allele.region.contig_region(), sequence: allele.sequence }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn genomic(contig: &str, begin: u32, end: u32) -> GenomicRegion {
        GenomicRegion::new(contig, begin, end).unwrap()
    }

    #[test]
    fn allele_from_parts() {
        let a = Allele::from_parts("chr1", 100, b"ACG".to_vec()).unwrap();
        assert_eq!(a.mapped_region().contig_name(), "chr1");
        assert_eq!(a.mapped_region().begin(), 100);
        assert_eq!(a.mapped_region().end(), 103);
        assert_eq!(a.sequence(), b"ACG");
        assert_eq!(a.sequence_size(), 3);
        assert!(!a.is_sequence_empty());
    }

    #[test]
    fn allele_from_parts_empty_sequence() {
        let a = Allele::from_parts("chrX", 50, vec![]).unwrap();
        assert_eq!(a.mapped_region().begin(), 50);
        assert_eq!(a.mapped_region().end(), 50);
        assert!(a.is_sequence_empty());
    }

    #[test]
    fn contig_allele_from_position() {
        let ca = ContigAllele::from_position(200, b"TTAG".to_vec()).unwrap();
        assert_eq!(ca.mapped_region().begin(), 200);
        assert_eq!(ca.mapped_region().end(), 204);
        assert_eq!(ca.sequence_size(), 4);
    }

    #[test]
    fn is_insertion_allele_sequence_longer_than_region() {
        let region = genomic("chr1", 100, 100);
        let a = Allele::new(region, b"ACG".to_vec());
        assert!(is_insertion_allele(&a));
        assert!(!is_deletion_allele(&a));
        assert!(is_indel_allele(&a));
    }

    #[test]
    fn is_deletion_allele_sequence_shorter_than_region() {
        let region = genomic("chr1", 100, 104);
        let a = Allele::new(region, vec![]);
        assert!(is_deletion_allele(&a));
        assert!(!is_insertion_allele(&a));
        assert!(is_indel_allele(&a));
    }

    #[test]
    fn substitution_allele_is_not_indel() {
        let region = genomic("chr1", 100, 101);
        let a = Allele::new(region, b"T".to_vec());
        assert!(!is_indel_allele(&a));
        assert!(!is_insertion_allele(&a));
        assert!(!is_deletion_allele(&a));
    }

    #[test]
    fn simple_insertion_empty_region_with_sequence() {
        let region = genomic("chr1", 100, 100);
        let a = Allele::new(region, b"ACG".to_vec());
        assert!(is_simple_insertion(&a));
    }

    #[test]
    fn not_simple_insertion_when_region_nonempty() {
        let region = genomic("chr1", 100, 101);
        let a = Allele::new(region, b"ACGT".to_vec());
        assert!(!is_simple_insertion(&a));
    }

    #[test]
    fn simple_deletion_nonempty_region_empty_sequence() {
        let region = genomic("chr1", 100, 105);
        let a = Allele::new(region, vec![]);
        assert!(is_simple_deletion(&a));
    }

    #[test]
    fn not_simple_deletion_when_sequence_nonempty() {
        let region = genomic("chr1", 100, 105);
        let a = Allele::new(region, b"A".to_vec());
        assert!(!is_simple_deletion(&a));
    }

    #[test]
    fn allele_ordering_by_region_then_sequence() {
        let a1 = Allele::from_parts("chr1", 100, b"A".to_vec()).unwrap();
        let a2 = Allele::from_parts("chr1", 100, b"T".to_vec()).unwrap();
        let a3 = Allele::from_parts("chr1", 200, b"A".to_vec()).unwrap();
        assert!(a1 < a2);
        assert!(a1 < a3);
    }

    #[test]
    fn allele_ordering_different_contigs_is_none() {
        let a1 = Allele::from_parts("chr1", 100, b"A".to_vec()).unwrap();
        let a2 = Allele::from_parts("chr2", 100, b"A".to_vec()).unwrap();
        assert_eq!(a1.partial_cmp(&a2), None);
    }

    #[test]
    fn contig_allele_ordering() {
        let ca1 = ContigAllele::from_position(100, b"A".to_vec()).unwrap();
        let ca2 = ContigAllele::from_position(200, b"A".to_vec()).unwrap();
        assert!(ca1 < ca2);
    }

    #[test]
    fn allele_equality() {
        let a1 = Allele::from_parts("chr1", 100, b"ACG".to_vec()).unwrap();
        let a2 = Allele::from_parts("chr1", 100, b"ACG".to_vec()).unwrap();
        let a3 = Allele::from_parts("chr1", 100, b"ACT".to_vec()).unwrap();
        assert_eq!(a1, a2);
        assert_ne!(a1, a3);
    }

    #[test]
    fn demote_strips_contig_name() {
        let a = Allele::from_parts("chr1", 100, b"ACG".to_vec()).unwrap();
        let ca = demote(a);
        assert_eq!(ca.mapped_region().begin(), 100);
        assert_eq!(ca.mapped_region().end(), 103);
        assert_eq!(ca.sequence(), b"ACG");
    }

    #[test]
    fn allele_display() {
        let a = Allele::from_parts("chr1", 10, b"AT".to_vec()).unwrap();
        let s = a.to_string();
        assert!(s.contains("AT"));
        assert!(s.contains("chr1"));
    }
}
