// Converted from C++ to Rust

use std::collections::HashMap;
use crate::io::variant::vcf_record::{VcfRecord, VcfRecordBuilder};
use crate::core::types::variant::Variant;
use crate::core::types::genotype::Genotype;
use crate::core::types::haplotype::Haplotype;
use crate::basics::genomic_region::GenomicRegion;

pub struct VcfRecordFactory;

impl VcfRecordFactory {
    pub fn make_record(
        variant: &Variant,
        sample_genotypes: &HashMap<String, (Genotype, f64)>,
        qual: Option<f32>,
    ) -> Result<VcfRecord, String> {
        let region = variant.mapped_region();
        let ref_seq = String::from_utf8_lossy(variant.ref_sequence()).to_string();
        let alt_seq = String::from_utf8_lossy(variant.alt_sequence()).to_string();

        let filters = if qual.unwrap_or(0.0) < 2.0 {
            vec!["LowQual".to_string()]
        } else {
            vec!["PASS".to_string()]
        };

        let info = HashMap::new();
        let mut record = VcfRecord::new(
            region.clone(),
            ".".to_string(),
            ref_seq,
            vec![alt_seq],
            qual,
            filters,
            info,
        );

        let format_keys = vec!["GT".to_string(), "GQ".to_string()];
        let mut samples = HashMap::new();

        for (sample, (genotype, posterior)) in sample_genotypes {
            let gt = genotype_to_gt_string(genotype, variant);
            let gq = (posterior * 100.0) as i64;
            let mut sample_data = HashMap::new();
            sample_data.insert("GT".to_string(), vec![gt]);
            sample_data.insert("GQ".to_string(), vec![gq.to_string()]);
            samples.insert(sample.clone(), sample_data);
        }

        Ok(record.with_genotypes(format_keys, samples))
    }
}

fn genotype_to_gt_string(genotype: &Genotype, variant: &Variant) -> String {
    let ploidy = genotype.ploidy();
    if ploidy == 0 { return ".".to_string(); }

    let allele_indices: Vec<String> = (0..ploidy).map(|i| {
        let h = &genotype.haplotypes()[i];
        let contains_alt = h.sequence().windows(variant.alt_sequence().len())
            .any(|w| w == variant.alt_sequence().as_slice());
        if contains_alt { "1".to_string() } else { "0".to_string() }
    }).collect();

    allele_indices.join("/")
}
