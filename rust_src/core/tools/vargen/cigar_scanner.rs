// Converted from C++ to Rust

use crate::core::types::variant::Variant;
use crate::basics::genomic_region::GenomicRegion;
use crate::basics::aligned_read::AlignedRead;
use crate::basics::cigar_string::{CigarFlag, CigarOperation};
use crate::io::reference::reference_genome::ReferenceGenome;
use super::variant_generator::VariantGenerator;

pub struct CigarScanner {
    min_base_quality: u8,
}

impl CigarScanner {
    pub fn new(min_base_quality: u8) -> Self {
        CigarScanner { min_base_quality }
    }
}

impl VariantGenerator for CigarScanner {
    fn name(&self) -> &str { "CigarScanner" }

    fn generate(&self, reads: &[AlignedRead], region: &GenomicRegion, reference: &ReferenceGenome) -> Vec<Variant> {
        let mut variants = Vec::new();

        let ref_seq = match reference.fetch_sequence(region) {
            Ok(s) => s,
            Err(_) => return variants,
        };

        for read in reads {
            if !read.mapped_region().overlaps(region) { continue; }
            let read_variants = self.scan_read(read, &ref_seq, region);
            variants.extend(read_variants);
        }

        variants.sort();
        variants.dedup();
        variants
    }
}

impl CigarScanner {
    fn scan_read(&self, read: &AlignedRead, ref_seq: &[u8], region: &GenomicRegion) -> Vec<Variant> {
        let mut variants = Vec::new();
        let cigar = read.cigar();
        let read_seq = read.sequence();
        let qualities = read.qualities();

        let mut ref_pos = read.mapped_region().begin();
        let mut read_pos = 0usize;
        let region_begin = region.begin();

        for op in cigar.iter() {
            match op.flag() {
                CigarFlag::AlignmentMatch | CigarFlag::SequenceMatch | CigarFlag::Substitution => {
                    for i in 0..op.size() as usize {
                        if ref_pos < region_begin {
                            ref_pos += 1;
                            read_pos += 1;
                            continue;
                        }
                        let ref_idx = (ref_pos - region_begin) as usize;
                        if ref_idx >= ref_seq.len() || read_pos >= read_seq.len() { break; }

                        let qual = qualities.get(read_pos).copied().unwrap_or(0);
                        if qual >= self.min_base_quality && ref_seq[ref_idx] != read_seq[read_pos] {
                            if let Ok(v) = Variant::from_parts(
                                region.contig_name(),
                                ref_pos,
                                vec![ref_seq[ref_idx]],
                                vec![read_seq[read_pos]],
                            ) {
                                variants.push(v);
                            }
                        }
                        ref_pos += 1;
                        read_pos += 1;
                    }
                }
                CigarFlag::Insertion => {
                    if ref_pos >= region_begin {
                        let end = (read_pos + op.size() as usize).min(read_seq.len());
                        let alt_seq: Vec<u8> = read_seq[read_pos..end].to_vec();
                        if let Ok(v) = Variant::from_parts(
                            region.contig_name(), ref_pos, vec![], alt_seq,
                        ) {
                            variants.push(v);
                        }
                    }
                    read_pos += op.size() as usize;
                }
                CigarFlag::Deletion => {
                    if ref_pos >= region_begin {
                        let ref_idx = (ref_pos - region_begin) as usize;
                        let del_len = op.size() as usize;
                        if ref_idx + del_len <= ref_seq.len() {
                            let del_seq = ref_seq[ref_idx..ref_idx + del_len].to_vec();
                            if let Ok(v) = Variant::from_parts(
                                region.contig_name(), ref_pos, del_seq, vec![],
                            ) {
                                variants.push(v);
                            }
                        }
                    }
                    ref_pos += op.size();
                }
                CigarFlag::SoftClipped => {
                    read_pos += op.size() as usize;
                }
                CigarFlag::HardClipped | CigarFlag::Padding => {}
                _ => {
                    ref_pos += op.size();
                    read_pos += op.size() as usize;
                }
            }
        }

        variants
    }
}
