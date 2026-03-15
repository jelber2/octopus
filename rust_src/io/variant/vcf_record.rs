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

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn region(contig: &str, begin: u32, end: u32) -> GenomicRegion {
        GenomicRegion::new(contig, begin, end).unwrap()
    }

    fn simple_record() -> VcfRecord {
        VcfRecord::new(
            region("chr1", 99, 100),
            ".".to_string(),
            "A".to_string(),
            vec!["T".to_string()],
            Some(50.0),
            vec!["PASS".to_string()],
            HashMap::new(),
        )
    }

    #[test]
    fn chrom_and_pos_1_based() {
        let rec = simple_record();
        assert_eq!(rec.chrom(), "chr1");
        assert_eq!(rec.pos(), 100);
    }

    #[test]
    fn ref_allele_and_alt() {
        let rec = simple_record();
        assert_eq!(rec.ref_allele(), "A");
        assert_eq!(rec.num_alt(), 1);
        assert_eq!(rec.alt()[0], "T");
    }

    #[test]
    fn qual_is_some() {
        let rec = simple_record();
        assert!((rec.qual().unwrap() - 50.0_f32).abs() < 1e-6);
    }

    #[test]
    fn qual_is_none_when_not_set() {
        let rec = VcfRecord::new(
            region("chr1", 99, 100),
            ".".to_string(),
            "A".to_string(),
            vec!["T".to_string()],
            None,
            vec![],
            HashMap::new(),
        );
        assert!(rec.qual().is_none());
    }

    #[test]
    fn is_pass_with_pass_filter() {
        let rec = simple_record();
        assert!(rec.is_pass());
    }

    #[test]
    fn is_pass_with_empty_filters() {
        let rec = VcfRecord::new(
            region("chr1", 0, 1),
            ".".to_string(),
            "A".to_string(),
            vec!["T".to_string()],
            None,
            vec![],
            HashMap::new(),
        );
        assert!(rec.is_pass());
    }

    #[test]
    fn is_not_pass_with_fail_filter() {
        let rec = VcfRecord::new(
            region("chr1", 0, 1),
            ".".to_string(),
            "A".to_string(),
            vec!["T".to_string()],
            None,
            vec!["LowQual".to_string()],
            HashMap::new(),
        );
        assert!(!rec.is_pass());
        assert!(rec.has_filter("LowQual"));
        assert!(!rec.has_filter("PASS"));
    }

    #[test]
    fn has_filter_checks_membership() {
        let mut rec = simple_record();
        assert!(rec.has_filter("PASS"));
        assert!(!rec.has_filter("LowQual"));
    }

    #[test]
    fn info_key_access() {
        let mut info = HashMap::new();
        info.insert("DP".to_string(), vec!["100".to_string()]);
        let rec = VcfRecord::new(
            region("chr1", 0, 1),
            ".".to_string(),
            "A".to_string(),
            vec!["T".to_string()],
            None,
            vec!["PASS".to_string()],
            info,
        );
        assert!(rec.has_info("DP"));
        assert!(!rec.has_info("AF"));
        let dp = rec.info_value("DP").unwrap();
        assert_eq!(dp[0], "100");
    }

    #[test]
    fn with_genotypes_adds_format_and_samples() {
        let mut samples: HashMap<String, HashMap<String, Vec<String>>> = HashMap::new();
        let mut sample_data = HashMap::new();
        sample_data.insert("GT".to_string(), vec!["0/1".to_string()]);
        sample_data.insert("GQ".to_string(), vec!["30".to_string()]);
        samples.insert("SAMPLE1".to_string(), sample_data);

        let rec = simple_record().with_genotypes(
            vec!["GT".to_string(), "GQ".to_string()],
            samples,
        );

        assert!(rec.has_format("GT"));
        assert!(rec.has_format("GQ"));
        assert!(!rec.has_format("DP"));
        assert!(rec.has_sample("SAMPLE1"));
        assert!(!rec.has_sample("SAMPLE2"));
    }

    #[test]
    fn genotype_heterozygous_unphased() {
        let mut samples: HashMap<String, HashMap<String, Vec<String>>> = HashMap::new();
        let mut sd = HashMap::new();
        sd.insert("GT".to_string(), vec!["0/1".to_string()]);
        samples.insert("S1".to_string(), sd);

        let rec = simple_record().with_genotypes(vec!["GT".to_string()], samples);
        let gt = rec.genotype("S1").unwrap();
        assert_eq!(gt, vec![0, 1]);
    }

    #[test]
    fn genotype_homozygous_alt() {
        let mut samples: HashMap<String, HashMap<String, Vec<String>>> = HashMap::new();
        let mut sd = HashMap::new();
        sd.insert("GT".to_string(), vec!["1/1".to_string()]);
        samples.insert("S1".to_string(), sd);

        let rec = simple_record().with_genotypes(vec!["GT".to_string()], samples);
        let gt = rec.genotype("S1").unwrap();
        assert_eq!(gt, vec![1, 1]);
    }

    #[test]
    fn genotype_phased() {
        let mut samples: HashMap<String, HashMap<String, Vec<String>>> = HashMap::new();
        let mut sd = HashMap::new();
        sd.insert("GT".to_string(), vec!["0|1".to_string()]);
        samples.insert("S1".to_string(), sd);

        let rec = simple_record().with_genotypes(vec!["GT".to_string()], samples);
        let gt = rec.genotype("S1").unwrap();
        assert_eq!(gt, vec![0, 1]);
    }

    #[test]
    fn genotype_missing_allele() {
        let mut samples: HashMap<String, HashMap<String, Vec<String>>> = HashMap::new();
        let mut sd = HashMap::new();
        sd.insert("GT".to_string(), vec!["./1".to_string()]);
        samples.insert("S1".to_string(), sd);

        let rec = simple_record().with_genotypes(vec!["GT".to_string()], samples);
        let gt = rec.genotype("S1").unwrap();
        assert_eq!(gt, vec![-1, 1]);
    }

    #[test]
    fn genotype_missing_sample_returns_none() {
        let rec = simple_record();
        assert!(rec.genotype("NONEXISTENT").is_none());
    }

    #[test]
    fn get_sample_value_retrieves_format_field() {
        let mut samples: HashMap<String, HashMap<String, Vec<String>>> = HashMap::new();
        let mut sd = HashMap::new();
        sd.insert("DP".to_string(), vec!["45".to_string()]);
        samples.insert("S1".to_string(), sd);

        let rec = simple_record().with_genotypes(vec!["DP".to_string()], samples);
        let dp = rec.get_sample_value("S1", "DP").unwrap();
        assert_eq!(dp[0], "45");
    }

    #[test]
    fn display_shows_chrom_pos_ref_alt() {
        let rec = simple_record();
        let s = rec.to_string();
        assert!(s.contains("chr1"));
        assert!(s.contains("100"));
        assert!(s.contains("A"));
        assert!(s.contains("T"));
    }

    #[test]
    fn multiple_alt_alleles() {
        let rec = VcfRecord::new(
            region("chr1", 99, 100),
            "rs123".to_string(),
            "A".to_string(),
            vec!["T".to_string(), "C".to_string(), "G".to_string()],
            Some(99.0),
            vec!["PASS".to_string()],
            HashMap::new(),
        );
        assert_eq!(rec.num_alt(), 3);
        assert_eq!(rec.id(), "rs123");
    }

    #[test]
    fn builder_constructs_record() {
        let rec = VcfRecordBuilder::new()
            .region(region("chr2", 999, 1000))
            .id("rs456")
            .ref_allele("G")
            .add_alt("A")
            .qual(75.0)
            .add_filter("PASS")
            .add_info("AF", vec!["0.5".to_string()])
            .build()
            .unwrap();

        assert_eq!(rec.chrom(), "chr2");
        assert_eq!(rec.pos(), 1000);
        assert_eq!(rec.id(), "rs456");
        assert_eq!(rec.ref_allele(), "G");
        assert_eq!(rec.alt()[0], "A");
        assert!((rec.qual().unwrap() - 75.0_f32).abs() < 1e-6);
        assert!(rec.has_filter("PASS"));
        assert!(rec.has_info("AF"));
    }

    #[test]
    fn builder_without_region_errors() {
        let result = VcfRecordBuilder::new()
            .ref_allele("A")
            .add_alt("T")
            .build();
        assert!(result.is_err());
    }
}
