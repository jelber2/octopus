// Converted from C++ to Rust

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::collections::HashMap;
use super::vcf_record::{VcfRecord, VcfRecordBuilder, NucleotideSequence, KeyType, ValueType};
use super::vcf_header::VcfHeader;
use crate::basics::genomic_region::GenomicRegion;

pub struct VcfReader {
    header: VcfHeader,
    records: Vec<VcfRecord>,
}

impl VcfReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, String> {
        let file = File::open(path.as_ref()).map_err(|e| e.to_string())?;
        let reader = BufReader::new(file);
        let mut header = VcfHeader::new("VCFv4.3");
        let mut records = Vec::new();
        let mut samples: Vec<String> = Vec::new();
        let mut in_header = true;

        for line in reader.lines() {
            let line = line.map_err(|e| e.to_string())?;
            if line.starts_with("##") {
                parse_meta_line(&line, &mut header);
            } else if line.starts_with('#') {
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() > 9 {
                    samples = fields[9..].iter().map(|s| s.to_string()).collect();
                    for s in &samples {
                        header.add_sample(s.clone());
                    }
                }
                in_header = false;
            } else if !line.is_empty() {
                if let Ok(record) = parse_record_line(&line, &samples) {
                    records.push(record);
                }
            }
        }

        Ok(VcfReader { header, records })
    }

    pub fn header(&self) -> &VcfHeader { &self.header }
    pub fn records(&self) -> &[VcfRecord] { &self.records }

    pub fn fetch(&self, region: &GenomicRegion) -> Vec<&VcfRecord> {
        self.records.iter().filter(|r| {
            r.chrom() == region.contig_name()
                && r.mapped_region().begin() < region.end()
                && r.mapped_region().end() > region.begin()
        }).collect()
    }
}

fn parse_meta_line(line: &str, header: &mut VcfHeader) {
    if let Some(rest) = line.strip_prefix("##contig=<") {
        let id = extract_tag(rest, "ID").unwrap_or_default();
        let length = extract_tag(rest, "length").and_then(|l| l.parse().ok());
        header.add_contig(id, length);
    } else if let Some(rest) = line.strip_prefix("##FILTER=<") {
        let id = extract_tag(rest, "ID").unwrap_or_default();
        let desc = extract_tag(rest, "Description").unwrap_or_default();
        header.add_filter(id, desc);
    } else if let Some(rest) = line.strip_prefix("##INFO=<") {
        let id = extract_tag(rest, "ID").unwrap_or_default();
        let num = extract_tag(rest, "Number").unwrap_or(".".to_string());
        let type_ = extract_tag(rest, "Type").unwrap_or("String".to_string());
        let desc = extract_tag(rest, "Description").unwrap_or_default();
        header.add_info(id, num, type_, desc);
    } else if let Some(rest) = line.strip_prefix("##FORMAT=<") {
        let id = extract_tag(rest, "ID").unwrap_or_default();
        let num = extract_tag(rest, "Number").unwrap_or(".".to_string());
        let type_ = extract_tag(rest, "Type").unwrap_or("String".to_string());
        let desc = extract_tag(rest, "Description").unwrap_or_default();
        header.add_format(id, num, type_, desc);
    }
}

fn extract_tag(s: &str, tag: &str) -> Option<String> {
    let search = format!("{}=", tag);
    let pos = s.find(&search)?;
    let rest = &s[pos + search.len()..];
    if rest.starts_with('"') {
        let end = rest[1..].find('"')?;
        Some(rest[1..end + 1].to_string())
    } else {
        let end = rest.find(|c| c == ',' || c == '>').unwrap_or(rest.len());
        Some(rest[..end].to_string())
    }
}

fn parse_record_line(line: &str, samples: &[String]) -> Result<VcfRecord, String> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 8 {
        return Err(format!("VCF record has too few fields: {}", line));
    }

    let chrom = fields[0].to_string();
    let pos: u32 = fields[1].parse().map_err(|e: std::num::ParseIntError| e.to_string())?;
    let id = fields[2].to_string();
    let ref_allele = fields[3].to_string();
    let alt_alleles: Vec<String> = fields[4].split(',').map(|s| s.to_string()).collect();
    let qual: Option<f32> = if fields[5] == "." { None } else { fields[5].parse().ok() };
    let filters: Vec<String> = if fields[6] == "." || fields[6] == "PASS" {
        vec![fields[6].to_string()]
    } else {
        fields[6].split(';').map(|s| s.to_string()).collect()
    };

    let mut info = HashMap::new();
    if fields[7] != "." {
        for kv in fields[7].split(';') {
            let parts: Vec<&str> = kv.splitn(2, '=').collect();
            let key = parts[0].to_string();
            let values = if parts.len() > 1 {
                parts[1].split(',').map(|s| s.to_string()).collect()
            } else {
                vec!["true".to_string()]
            };
            info.insert(key, values);
        }
    }

    let begin = pos.saturating_sub(1);
    let end = begin + ref_allele.len() as u32;
    let region = GenomicRegion::new(chrom, begin, end).map_err(|e| e.to_string())?;

    let mut record = VcfRecord::new(region, id, ref_allele, alt_alleles, qual, filters, info);

    if fields.len() > 8 && !samples.is_empty() {
        let format_keys: Vec<String> = fields[8].split(':').map(|s| s.to_string()).collect();
        let mut sample_map = HashMap::new();

        for (i, sample) in samples.iter().enumerate() {
            if let Some(sample_field) = fields.get(9 + i) {
                let values: Vec<&str> = sample_field.split(':').collect();
                let mut sample_data = HashMap::new();
                for (j, key) in format_keys.iter().enumerate() {
                    let vals = values.get(j).unwrap_or(&".");
                    let parsed: Vec<String> = vals.split(',').map(|s| s.to_string()).collect();
                    sample_data.insert(key.clone(), parsed);
                }
                sample_map.insert(sample.clone(), sample_data);
            }
        }
        record = record.with_genotypes(format_keys, sample_map);
    }

    Ok(record)
}
