// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Converted from C++ to Rust

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use crate::basics::genomic_region::GenomicRegion;
use crate::basics::contig_region::{ContigRegion, Position, Size};

pub type ContigName = String;
pub type GeneticSequence = Vec<u8>;

pub trait ReferenceReader: Send + Sync {
    fn name(&self) -> &str;
    fn contig_names(&self) -> Vec<ContigName>;
    fn contig_size(&self, contig: &str) -> Option<Size>;
    fn fetch_sequence(&self, region: &GenomicRegion) -> Result<GeneticSequence, String>;
}

pub struct ReferenceGenome {
    impl_: Box<dyn ReferenceReader>,
    name: String,
    contig_sizes: HashMap<ContigName, Size>,
    ordered_contigs: Vec<ContigName>,
}

impl ReferenceGenome {
    pub fn new(impl_: Box<dyn ReferenceReader>) -> Self {
        let name = impl_.name().to_string();
        let contig_names = impl_.contig_names();
        let mut contig_sizes = HashMap::new();
        for cn in &contig_names {
            if let Some(size) = impl_.contig_size(cn) {
                contig_sizes.insert(cn.clone(), size);
            }
        }
        ReferenceGenome {
            impl_,
            name,
            contig_sizes,
            ordered_contigs: contig_names,
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn has_contig(&self, contig: &str) -> bool {
        self.contig_sizes.contains_key(contig)
    }

    pub fn num_contigs(&self) -> usize {
        self.ordered_contigs.len()
    }

    pub fn contig_names(&self) -> &[ContigName] {
        &self.ordered_contigs
    }

    pub fn contig_size(&self, contig: &str) -> Option<Size> {
        self.contig_sizes.get(contig).copied()
    }

    pub fn contig_region(&self, contig: &str) -> Option<GenomicRegion> {
        let size = self.contig_size(contig)?;
        GenomicRegion::new(contig, 0, size).ok()
    }

    pub fn contains(&self, region: &GenomicRegion) -> bool {
        if let Some(size) = self.contig_size(region.contig_name()) {
            region.end() <= size
        } else {
            false
        }
    }

    pub fn fetch_sequence(&self, region: &GenomicRegion) -> Result<GeneticSequence, String> {
        self.impl_.fetch_sequence(region)
    }
}

pub fn get_all_contig_regions(reference: &ReferenceGenome) -> Vec<GenomicRegion> {
    reference.contig_names().iter()
        .filter_map(|cn| reference.contig_region(cn))
        .collect()
}

pub fn calculate_genome_size(reference: &ReferenceGenome) -> Position {
    reference.contig_names().iter()
        .filter_map(|cn| reference.contig_size(cn))
        .sum()
}
