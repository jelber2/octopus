// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Converted from C++ to Rust

use std::collections::HashMap;
use std::fmt;
use crate::basics::genomic_region::GenomicRegion;

pub type NucleotideSequence = String;
pub type QualityType = f32;
pub type SampleName = String;
pub type KeyType = String;
pub type ValueType = String;

#[derive(Debug, Clone)]
pub struct VcfRecord {
    region: GenomicRegion,
    id: String,
    ref_allele: NucleotideSequence,
    alt_alleles: Vec<NucleotideSequence>,
    qual: Option<QualityType>,
    filters: Vec<KeyType>,
    info: HashMap<KeyType, Vec<ValueType>>,
    format: Vec<KeyType>,
    samples: HashMap<SampleName, HashMap<KeyType, Vec<ValueType>>>,
}

impl VcfRecord {
    pub fn new(
        region: GenomicRegion,
        id: String,
        ref_allele: NucleotideSequence,
        alt_alleles: Vec<NucleotideSequence>,
        qual: Option<QualityType>,
        filters: Vec<KeyType>,
        info: HashMap<KeyType, Vec<ValueType>>,
    ) -> Self {
        VcfRecord { region, id, ref_allele, alt_alleles, qual, filters, info, format: vec![], samples: HashMap::new() }
    }

    pub fn with_genotypes(mut self, format: Vec<KeyType>, samples: HashMap<SampleName, HashMap<KeyType, Vec<ValueType>>>) -> Self {
        self.format = format;
        self.samples = samples;
        self
    }

    pub fn mapped_region(&self) -> &GenomicRegion { &self.region }
    pub fn chrom(&self) -> &str { self.region.contig_name() }
    pub fn pos(&self) -> u32 { self.region.begin() + 1 }
    pub fn id(&self) -> &str { &self.id }
    pub fn ref_allele(&self) -> &NucleotideSequence { &self.ref_allele }
    pub fn num_alt(&self) -> usize { self.alt_alleles.len() }
    pub fn alt(&self) -> &[NucleotideSequence] { &self.alt_alleles }
    pub fn qual(&self) -> Option<QualityType> { self.qual }

    pub fn has_filter(&self, filter: &str) -> bool {
        self.filters.iter().any(|f| f == filter)
    }
    pub fn filters(&self) -> &[KeyType] { &self.filters }
    pub fn is_pass(&self) -> bool {
        self.filters.is_empty() || (self.filters.len() == 1 && self.filters[0] == "PASS")
    }

    pub fn has_info(&self, key: &str) -> bool {
        self.info.contains_key(key)
    }
    pub fn info_keys(&self) -> Vec<&KeyType> {
        self.info.keys().collect()
    }
    pub fn info_value(&self, key: &str) -> Option<&Vec<ValueType>> {
        self.info.get(key)
    }

    pub fn has_format(&self, key: &str) -> bool {
        self.format.iter().any(|f| f == key)
    }
    pub fn format(&self) -> &[KeyType] { &self.format }

    pub fn has_sample(&self, sample: &str) -> bool {
        self.samples.contains_key(sample)
    }
    pub fn samples(&self) -> Vec<&SampleName> {
        self.samples.keys().collect()
    }
    pub fn get_sample_value(&self, sample: &str, key: &str) -> Option<&Vec<ValueType>> {
        self.samples.get(sample)?.get(key)
    }

    pub fn genotype(&self, sample: &str) -> Option<Vec<i8>> {
        let gt_values = self.get_sample_value(sample, "GT")?;
        let gt = gt_values.first()?;
        let phased = gt.contains('|');
        let sep = if phased { '|' } else { '/' };
        Some(gt.split(sep).map(|a| {
            if a == "." { -1i8 } else { a.parse().unwrap_or(-1) }
        }).collect())
    }
}

impl fmt::Display for VcfRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}\t{}\t{}\t{}\t{}",
            self.chrom(), self.pos(), self.id(),
            self.ref_allele, self.alt_alleles.join(","))
    }
}

pub struct VcfRecordBuilder {
    region: Option<GenomicRegion>,
    id: String,
    ref_allele: NucleotideSequence,
    alt_alleles: Vec<NucleotideSequence>,
    qual: Option<QualityType>,
    filters: Vec<KeyType>,
    info: HashMap<KeyType, Vec<ValueType>>,
    format: Vec<KeyType>,
    samples: HashMap<SampleName, HashMap<KeyType, Vec<ValueType>>>,
}

impl VcfRecordBuilder {
    pub fn new() -> Self {
        VcfRecordBuilder {
            region: None, id: ".".to_string(), ref_allele: String::new(),
            alt_alleles: vec![], qual: None, filters: vec![], info: HashMap::new(),
            format: vec![], samples: HashMap::new(),
        }
    }

    pub fn region(mut self, region: GenomicRegion) -> Self { self.region = Some(region); self }
    pub fn id(mut self, id: impl Into<String>) -> Self { self.id = id.into(); self }
    pub fn ref_allele(mut self, r: impl Into<String>) -> Self { self.ref_allele = r.into(); self }
    pub fn add_alt(mut self, alt: impl Into<String>) -> Self { self.alt_alleles.push(alt.into()); self }
    pub fn qual(mut self, q: f32) -> Self { self.qual = Some(q); self }
    pub fn add_filter(mut self, f: impl Into<String>) -> Self { self.filters.push(f.into()); self }
    pub fn add_info(mut self, k: impl Into<String>, v: Vec<String>) -> Self { self.info.insert(k.into(), v); self }

    pub fn build(self) -> Result<VcfRecord, String> {
        let region = self.region.ok_or("region not set")?;
        Ok(VcfRecord {
            region, id: self.id, ref_allele: self.ref_allele,
            alt_alleles: self.alt_alleles, qual: self.qual, filters: self.filters,
            info: self.info, format: self.format, samples: self.samples,
        })
    }
}

impl Default for VcfRecordBuilder {
    fn default() -> Self { Self::new() }
}
