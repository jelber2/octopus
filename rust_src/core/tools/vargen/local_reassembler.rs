// Converted from C++ to Rust

use crate::core::types::variant::Variant;
use crate::basics::genomic_region::GenomicRegion;
use crate::basics::aligned_read::AlignedRead;
use crate::io::reference::reference_genome::ReferenceGenome;
use super::variant_generator::VariantGenerator;
use std::collections::HashMap;

pub struct LocalReassembler {
    kmer_size: usize,
    min_coverage: usize,
}

impl LocalReassembler {
    pub fn new(kmer_size: usize, min_coverage: usize) -> Self {
        LocalReassembler { kmer_size, min_coverage }
    }

    fn build_de_bruijn_graph(&self, reads: &[AlignedRead]) -> HashMap<Vec<u8>, Vec<Vec<u8>>> {
        let mut graph: HashMap<Vec<u8>, Vec<Vec<u8>>> = HashMap::new();
        let k = self.kmer_size;

        for read in reads {
            let seq = read.sequence();
            if seq.len() < k { continue; }
            for i in 0..seq.len() - k {
                let left = seq[i..i + k - 1].to_vec();
                let right = seq[i + 1..i + k].to_vec();
                graph.entry(left).or_default().push(right);
            }
        }
        graph
    }
}

impl VariantGenerator for LocalReassembler {
    fn name(&self) -> &str { "LocalReassembler" }

    fn generate(&self, reads: &[AlignedRead], region: &GenomicRegion, reference: &ReferenceGenome) -> Vec<Variant> {
        Vec::new()
    }
}
