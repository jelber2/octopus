// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Converted from C++ to Rust

use std::fmt;
use std::hash::{Hash, Hasher};
use super::haplotype::Haplotype;
use super::allele::{Allele, NucleotideSequence};
use crate::basics::genomic_region::GenomicRegion;

#[derive(Debug, Clone)]
pub struct Genotype {
    haplotypes: Vec<Haplotype>,
}

impl Genotype {
    pub fn new(haplotypes: Vec<Haplotype>) -> Self {
        Genotype { haplotypes }
    }

    pub fn ploidy(&self) -> usize { self.haplotypes.len() }
    pub fn haplotypes(&self) -> &[Haplotype] { &self.haplotypes }
    pub fn is_empty(&self) -> bool { self.haplotypes.is_empty() }

    pub fn is_homozygous(&self) -> bool {
        if self.haplotypes.is_empty() { return true; }
        let first = &self.haplotypes[0];
        self.haplotypes.iter().all(|h| h == first)
    }

    pub fn is_heterozygous(&self) -> bool {
        !self.is_homozygous()
    }

    pub fn count(&self, haplotype: &Haplotype) -> usize {
        self.haplotypes.iter().filter(|h| *h == haplotype).count()
    }

    pub fn unique_haplotypes(&self) -> Vec<&Haplotype> {
        let mut seen = Vec::new();
        let mut result = Vec::new();
        for h in &self.haplotypes {
            if !seen.contains(&h.get_hash()) {
                seen.push(h.get_hash());
                result.push(h);
            }
        }
        result
    }
}

impl PartialEq for Genotype {
    fn eq(&self, other: &Self) -> bool {
        if self.ploidy() != other.ploidy() { return false; }
        let mut sorted_self: Vec<_> = self.haplotypes.iter().map(|h| h.get_hash()).collect();
        let mut sorted_other: Vec<_> = other.haplotypes.iter().map(|h| h.get_hash()).collect();
        sorted_self.sort();
        sorted_other.sort();
        sorted_self == sorted_other
    }
}

impl Eq for Genotype {}

impl Hash for Genotype {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let mut hashes: Vec<u64> = self.haplotypes.iter().map(|h| h.get_hash()).collect();
        hashes.sort();
        hashes.hash(state);
    }
}

impl fmt::Display for Genotype {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let seqs: Vec<String> = self.haplotypes.iter()
            .map(|h| String::from_utf8_lossy(h.sequence()).to_string())
            .collect();
        write!(f, "{}", seqs.join("/"))
    }
}

pub fn make_genotypes(haplotypes: &[Haplotype], ploidy: usize) -> Vec<Genotype> {
    if haplotypes.is_empty() || ploidy == 0 { return vec![]; }

    let n = haplotypes.len();
    let mut result = Vec::new();

    fn combinations(haplotypes: &[Haplotype], ploidy: usize, start: usize, current: &mut Vec<usize>, result: &mut Vec<Genotype>) {
        if current.len() == ploidy {
            let gt = Genotype::new(current.iter().map(|&i| haplotypes[i].clone()).collect());
            result.push(gt);
            return;
        }
        for i in start..haplotypes.len() {
            current.push(i);
            combinations(haplotypes, ploidy, i, current, result);
            current.pop();
        }
    }

    combinations(haplotypes, ploidy, 0, &mut vec![], &mut result);
    result
}
