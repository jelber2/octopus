// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Converted from C++ to Rust

use std::fmt;
use std::hash::{Hash, Hasher};
use crate::basics::genomic_region::GenomicRegion;
use super::allele::{Allele, NucleotideSequence};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant {
    reference: Allele,
    alternative: Allele,
}

#[derive(Debug, Clone)]
pub struct VariantError(pub String);

impl fmt::Display for VariantError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Variant error: {}", self.0)
    }
}

impl std::error::Error for VariantError {}

impl Variant {
    pub fn new(reference: Allele, alternative: Allele) -> Result<Self, VariantError> {
        if reference.region != alternative.region {
            return Err(VariantError("reference & alternative alleles must define the same region".to_string()));
        }
        Ok(Variant { reference, alternative })
    }

    pub fn from_parts(
        contig: impl Into<String>,
        begin: u32,
        ref_sequence: impl Into<NucleotideSequence>,
        alt_sequence: impl Into<NucleotideSequence>,
    ) -> Result<Self, String> {
        let contig = contig.into();
        let ref_seq = ref_sequence.into();
        let alt_seq = alt_sequence.into();
        let end = begin + ref_seq.len() as u32;

        let region = GenomicRegion::new(&contig, begin, end).map_err(|e| e.to_string())?;
        let reference = Allele::new(region.clone(), ref_seq);
        let alternative = Allele::new(region, alt_seq);
        Ok(Variant { reference, alternative })
    }

    pub fn mapped_region(&self) -> &GenomicRegion { &self.reference.region }
    pub fn ref_allele(&self) -> &Allele { &self.reference }
    pub fn alt_allele(&self) -> &Allele { &self.alternative }
    pub fn ref_sequence(&self) -> &NucleotideSequence { &self.reference.sequence }
    pub fn alt_sequence(&self) -> &NucleotideSequence { &self.alternative.sequence }
    pub fn ref_sequence_size(&self) -> usize { self.reference.sequence.len() }
    pub fn alt_sequence_size(&self) -> usize { self.alternative.sequence.len() }
}

impl PartialOrd for Variant {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Variant {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.reference.partial_cmp(&other.reference)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| self.alternative.partial_cmp(&other.alternative)
                .unwrap_or(std::cmp::Ordering::Equal))
    }
}

impl Hash for Variant {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.reference.hash(state);
        self.alternative.sequence.hash(state);
    }
}

impl fmt::Display for Variant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}/{}", self.mapped_region(),
            String::from_utf8_lossy(self.ref_sequence()),
            String::from_utf8_lossy(self.alt_sequence()))
    }
}

pub fn is_snv(v: &Variant) -> bool {
    v.ref_sequence_size() == 1 && v.alt_sequence_size() == 1
}

pub fn is_mnv(v: &Variant) -> bool {
    v.ref_sequence_size() > 1 && v.alt_sequence_size() == v.ref_sequence_size()
        && !is_snv(v)
}

pub fn is_insertion(v: &Variant) -> bool {
    v.alt_sequence_size() > v.ref_sequence_size()
}

pub fn is_deletion(v: &Variant) -> bool {
    v.alt_sequence_size() < v.ref_sequence_size()
}

pub fn is_indel(v: &Variant) -> bool {
    is_insertion(v) || is_deletion(v)
}

pub fn is_simple_insertion(v: &Variant) -> bool {
    v.ref_sequence_size() == 0 && v.alt_sequence_size() > 0
}

pub fn is_simple_deletion(v: &Variant) -> bool {
    v.alt_sequence_size() == 0 && v.ref_sequence_size() > 0
}

pub fn indel_size(v: &Variant) -> i64 {
    v.alt_sequence_size() as i64 - v.ref_sequence_size() as i64
}

pub fn is_parsimonious(v: &Variant) -> bool {
    if v.ref_sequence().is_empty() || v.alt_sequence().is_empty() {
        return false;
    }
    let min_len = v.ref_sequence_size().min(v.alt_sequence_size());
    if min_len == 0 { return false; }

    let front_match = v.ref_sequence().iter().zip(v.alt_sequence().iter())
        .take_while(|(r, a)| r == a).count();
    let back_match = v.ref_sequence().iter().rev().zip(v.alt_sequence().iter().rev())
        .take_while(|(r, a)| r == a).count();

    front_match == 0 && (min_len == 1 || back_match == 0)
}

pub fn remove_duplicates(variants: &mut Vec<Variant>) {
    variants.sort();
    variants.dedup();
}

pub fn split_mnv(variant: &Variant) -> Vec<Variant> {
    if !is_mnv(variant) { return vec![variant.clone()]; }
    let contig = variant.mapped_region().contig_name();
    let begin = variant.mapped_region().begin();
    variant.ref_sequence().iter().zip(variant.alt_sequence().iter()).enumerate()
        .filter(|(_, (r, a))| r != a)
        .filter_map(|(i, (r, a))| {
            Variant::from_parts(
                contig,
                begin + i as u32,
                vec![*r],
                vec![*a],
            ).ok()
        })
        .collect()
}
