// Converted from C++ to Rust

use std::fmt;
use super::genotype::Genotype;
use super::haplotype::Haplotype;

#[derive(Debug, Clone)]
pub struct CancerGenotype {
    germline: Genotype,
    somatic: Vec<Haplotype>,
}

impl CancerGenotype {
    pub fn new(germline: Genotype, somatic: Vec<Haplotype>) -> Self {
        CancerGenotype { germline, somatic }
    }

    pub fn germline(&self) -> &Genotype { &self.germline }
    pub fn somatic(&self) -> &[Haplotype] { &self.somatic }
    pub fn ploidy(&self) -> usize { self.germline.ploidy() + self.somatic.len() }

    pub fn is_somatic_empty(&self) -> bool { self.somatic.is_empty() }
}

impl fmt::Display for CancerGenotype {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Germline: {} | Somatic: {} haplotypes", self.germline, self.somatic.len())
    }
}
