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

#[cfg(test)]
mod tests {
    use super::*;

    fn snv(contig: &str, pos: u32, r: u8, a: u8) -> Variant {
        Variant::from_parts(contig, pos, vec![r], vec![a]).unwrap()
    }

    fn indel(contig: &str, pos: u32, r: &[u8], a: &[u8]) -> Variant {
        Variant::from_parts(contig, pos, r.to_vec(), a.to_vec()).unwrap()
    }

    #[test]
    fn from_parts_snv() {
        let v = snv("chr1", 100, b'A', b'T');
        assert_eq!(v.mapped_region().contig_name(), "chr1");
        assert_eq!(v.mapped_region().begin(), 100);
        assert_eq!(v.mapped_region().end(), 101);
        assert_eq!(v.ref_sequence(), b"A");
        assert_eq!(v.alt_sequence(), b"T");
    }

    #[test]
    fn from_parts_insertion() {
        let v = Variant::from_parts("chr1", 100, vec![], b"ACGT".to_vec()).unwrap();
        assert_eq!(v.ref_sequence_size(), 0);
        assert_eq!(v.alt_sequence_size(), 4);
    }

    #[test]
    fn from_parts_deletion() {
        let v = Variant::from_parts("chr1", 100, b"ACGT".to_vec(), vec![]).unwrap();
        assert_eq!(v.ref_sequence_size(), 4);
        assert_eq!(v.alt_sequence_size(), 0);
    }

    #[test]
    fn new_mismatched_regions_errors() {
        let region1 = GenomicRegion::new("chr1", 100, 101).unwrap();
        let region2 = GenomicRegion::new("chr1", 200, 201).unwrap();
        let ref_allele = Allele::new(region1, b"A".to_vec());
        let alt_allele = Allele::new(region2, b"T".to_vec());
        assert!(Variant::new(ref_allele, alt_allele).is_err());
    }

    #[test]
    fn is_snv_single_base_change() {
        let v = snv("chr1", 100, b'A', b'T');
        assert!(is_snv(&v));
        assert!(!is_mnv(&v));
        assert!(!is_indel(&v));
    }

    #[test]
    fn is_mnv_multi_base_same_length() {
        let v = indel("chr1", 100, b"ACG", b"TTT");
        assert!(is_mnv(&v));
        assert!(!is_snv(&v));
        assert!(!is_indel(&v));
    }

    #[test]
    fn insertion_has_longer_alt() {
        let v = Variant::from_parts("chr1", 100, vec![], b"ACG".to_vec()).unwrap();
        assert!(is_insertion(&v));
        assert!(!is_deletion(&v));
        assert!(is_indel(&v));
    }

    #[test]
    fn deletion_has_longer_ref() {
        let v = Variant::from_parts("chr1", 100, b"ACG".to_vec(), vec![]).unwrap();
        assert!(is_deletion(&v));
        assert!(!is_insertion(&v));
        assert!(is_indel(&v));
    }

    #[test]
    fn simple_insertion_has_zero_length_ref() {
        let v = Variant::from_parts("chr1", 100, vec![], b"ACG".to_vec()).unwrap();
        assert!(is_simple_insertion(&v));
        assert!(!is_simple_deletion(&v));
    }

    #[test]
    fn simple_deletion_has_zero_length_alt() {
        let v = Variant::from_parts("chr1", 100, b"ACG".to_vec(), vec![]).unwrap();
        assert!(is_simple_deletion(&v));
        assert!(!is_simple_insertion(&v));
    }

    #[test]
    fn indel_size_positive_for_insertion_negative_for_deletion() {
        let ins = Variant::from_parts("chr1", 100, b"A".to_vec(), b"ACGT".to_vec()).unwrap();
        assert_eq!(indel_size(&ins), 3);
        let del = Variant::from_parts("chr1", 100, b"ACGT".to_vec(), b"A".to_vec()).unwrap();
        assert_eq!(indel_size(&del), -3);
        let sub = snv("chr1", 100, b'A', b'T');
        assert_eq!(indel_size(&sub), 0);
    }

    #[test]
    fn is_parsimonious_no_padding_bases() {
        let v = snv("chr1", 100, b'A', b'T');
        assert!(is_parsimonious(&v));
    }

    #[test]
    fn is_not_parsimonious_shared_prefix() {
        let v = indel("chr1", 100, b"GAC", b"GTC");
        assert!(!is_parsimonious(&v));
    }

    #[test]
    fn ordering_by_position() {
        let a = snv("chr1", 100, b'A', b'T');
        let b = snv("chr1", 200, b'C', b'G');
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn ordering_same_position_by_sequence() {
        let a = Variant::from_parts("chr1", 100, b"A".to_vec(), b"C".to_vec()).unwrap();
        let b = Variant::from_parts("chr1", 100, b"A".to_vec(), b"T".to_vec()).unwrap();
        assert!(a < b);
    }

    #[test]
    fn equality() {
        let a = snv("chr1", 100, b'A', b'T');
        let b = snv("chr1", 100, b'A', b'T');
        let c = snv("chr1", 100, b'A', b'C');
        assert_eq!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn display_format() {
        let v = snv("chr1", 100, b'A', b'T');
        let s = v.to_string();
        assert!(s.contains("chr1"));
        assert!(s.contains("A/T"));
    }

    #[test]
    fn remove_duplicates_deduplicates_and_sorts() {
        let mut variants = vec![
            snv("chr1", 200, b'C', b'G'),
            snv("chr1", 100, b'A', b'T'),
            snv("chr1", 100, b'A', b'T'),
        ];
        remove_duplicates(&mut variants);
        assert_eq!(variants.len(), 2);
        assert_eq!(variants[0].mapped_region().begin(), 100);
        assert_eq!(variants[1].mapped_region().begin(), 200);
    }

    #[test]
    fn split_mnv_splits_into_snvs() {
        let v = indel("chr1", 100, b"ACG", b"ATT");
        let snvs = split_mnv(&v);
        assert_eq!(snvs.len(), 2);
        assert_eq!(snvs[0].ref_sequence(), b"C");
        assert_eq!(snvs[0].alt_sequence(), b"T");
        assert_eq!(snvs[1].ref_sequence(), b"G");
        assert_eq!(snvs[1].alt_sequence(), b"T");
    }

    #[test]
    fn split_non_mnv_returns_self() {
        let v = snv("chr1", 100, b'A', b'T');
        let result = split_mnv(&v);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], v);
    }

    #[test]
    fn ref_and_alt_accessors() {
        let v = snv("chr1", 50, b'G', b'C');
        assert_eq!(v.ref_sequence_size(), 1);
        assert_eq!(v.alt_sequence_size(), 1);
        assert_eq!(v.ref_allele().sequence(), b"G");
        assert_eq!(v.alt_allele().sequence(), b"C");
    }
}
