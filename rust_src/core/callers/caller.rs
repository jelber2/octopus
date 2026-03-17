// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Converted from C++ to Rust

use std::collections::HashMap;
use crate::basics::genomic_region::GenomicRegion;
use crate::io::variant::vcf_record::VcfRecord;
use crate::io::reference::reference_genome::ReferenceGenome;
use crate::basics::aligned_read::AlignedRead;

pub type SampleName = String;
pub type CallBuffer = Vec<VcfRecord>;

pub struct CallerEnvironment<'a> {
    pub reference: &'a ReferenceGenome,
    pub reads: HashMap<SampleName, Vec<AlignedRead>>,
    pub region: GenomicRegion,
}

pub trait Caller: Send {
    fn call_variants(&self, env: &CallerEnvironment) -> Result<CallBuffer, String>;
    fn name(&self) -> &str;
}

pub struct CallerOptions {
    pub min_variant_quality: f64,
    pub max_haplotypes: usize,
    pub min_read_depth: usize,
    pub ploidy: usize,
    /// Minimum base quality (Phred) a base must have to be used for candidate generation.
    pub min_base_quality: u8,
    /// Expected germline SNP heterozygosity (θ) for Hardy-Weinberg prior.
    pub snp_heterozygosity: f64,
    /// Expected germline indel heterozygosity (θ) for Hardy-Weinberg prior.
    pub indel_heterozygosity: f64,
}

impl Default for CallerOptions {
    fn default() -> Self {
        CallerOptions {
            min_variant_quality: 2.0,
            max_haplotypes: 128,
            min_read_depth: 1,
            ploidy: 2,
            min_base_quality: 20,
            snp_heterozygosity: 0.001,
            indel_heterozygosity: 0.0001,
        }
    }
}
