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
