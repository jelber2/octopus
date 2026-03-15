// Converted from C++ to Rust

use crate::core::types::variant::Variant;
use crate::basics::genomic_region::GenomicRegion;
use crate::basics::aligned_read::AlignedRead;
use crate::io::reference::reference_genome::ReferenceGenome;

pub trait VariantGenerator: Send + Sync {
    fn generate(&self, reads: &[AlignedRead], region: &GenomicRegion, reference: &ReferenceGenome) -> Vec<Variant>;
    fn name(&self) -> &str;
}
