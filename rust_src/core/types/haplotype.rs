// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Converted from C++ to Rust

use std::fmt;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
use crate::basics::genomic_region::GenomicRegion;
use crate::basics::contig_region::ContigRegion;
use crate::basics::cigar_string::CigarString;
use super::allele::{ContigAllele, Allele, NucleotideSequence};
use super::variant::Variant;

pub struct Haplotype {
    region: GenomicRegion,
    explicit_alleles: Vec<ContigAllele>,
    sequence: NucleotideSequence,
    cached_hash: u64,
}

impl Haplotype {
    pub fn new(region: GenomicRegion, sequence: NucleotideSequence) -> Self {
        let cached_hash = {
            let mut h = DefaultHasher::new();
            sequence.hash(&mut h);
            h.finish()
        };
        Haplotype {
            region,
            explicit_alleles: vec![],
            sequence,
            cached_hash,
        }
    }

    pub fn with_alleles(
        region: GenomicRegion,
        alleles: Vec<ContigAllele>,
        reference_sequence_fn: impl Fn(&GenomicRegion) -> NucleotideSequence,
    ) -> Self {
        let sequence = if alleles.is_empty() {
            reference_sequence_fn(&region)
        } else {
            build_sequence(&region, &alleles, reference_sequence_fn)
        };
        let cached_hash = {
            let mut h = DefaultHasher::new();
            sequence.hash(&mut h);
            h.finish()
        };
        Haplotype {
            region,
            explicit_alleles: alleles,
            sequence,
            cached_hash,
        }
    }

    pub fn mapped_region(&self) -> &GenomicRegion { &self.region }
    pub fn sequence(&self) -> &NucleotideSequence { &self.sequence }
    pub fn sequence_size(&self) -> usize { self.sequence.len() }

    pub fn contains_allele(&self, allele: &ContigAllele) -> bool {
        self.explicit_alleles.iter().any(|a| a == allele)
    }

    pub fn get_hash(&self) -> u64 { self.cached_hash }

    pub fn difference(&self, other: &Haplotype) -> Vec<Variant> {
        let mut result = Vec::new();
        if self.region != other.region { return result; }

        let s1 = &self.sequence;
        let s2 = &other.sequence;
        let begin = self.region.begin();
        let len = s1.len().min(s2.len());

        let mut i = 0;
        while i < len {
            if s1[i] != s2[i] {
                let start = i;
                while i < len && s1[i] != s2[i] { i += 1; }
                if let Ok(v) = Variant::from_parts(
                    self.region.contig_name(),
                    begin + start as u32,
                    s1[start..i].to_vec(),
                    s2[start..i].to_vec(),
                ) {
                    result.push(v);
                }
            } else {
                i += 1;
            }
        }
        result
    }
}

fn build_sequence(
    region: &GenomicRegion,
    alleles: &[ContigAllele],
    reference_fn: impl Fn(&GenomicRegion) -> NucleotideSequence,
) -> NucleotideSequence {
    let mut sequence = NucleotideSequence::new();
    let begin = region.begin();
    let end = region.end();

    let mut pos = begin;
    for allele in alleles {
        let allele_begin = allele.region.begin();
        if pos < allele_begin {
            if let Ok(ref_region) = GenomicRegion::new(region.contig_name(), pos, allele_begin) {
                sequence.extend(reference_fn(&ref_region));
            }
        }
        sequence.extend(&allele.sequence);
        pos = allele.region.end();
    }
    if pos < end {
        if let Ok(ref_region) = GenomicRegion::new(region.contig_name(), pos, end) {
            sequence.extend(reference_fn(&ref_region));
        }
    }
    sequence
}

impl PartialEq for Haplotype {
    fn eq(&self, other: &Self) -> bool {
        self.sequence == other.sequence && self.region == other.region
    }
}

impl Eq for Haplotype {}

impl Hash for Haplotype {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u64(self.cached_hash);
    }
}

impl fmt::Debug for Haplotype {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Haplotype({}, {})", self.region, String::from_utf8_lossy(&self.sequence))
    }
}

impl fmt::Display for Haplotype {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.region, String::from_utf8_lossy(&self.sequence))
    }
}

impl Clone for Haplotype {
    fn clone(&self) -> Self {
        Haplotype {
            region: self.region.clone(),
            explicit_alleles: self.explicit_alleles.clone(),
            sequence: self.sequence.clone(),
            cached_hash: self.cached_hash,
        }
    }
}

pub fn is_reference(haplotype: &Haplotype) -> bool {
    haplotype.explicit_alleles.is_empty()
}
