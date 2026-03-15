// Converted from C++ to Rust

use std::any::Any;
use super::facet::Facet;
use crate::core::types::allele::Allele;

pub struct AllelesFacet {
    alleles: Vec<Allele>,
}

impl AllelesFacet {
    pub fn new(alleles: Vec<Allele>) -> Self { AllelesFacet { alleles } }
    pub fn alleles(&self) -> &[Allele] { &self.alleles }
}

impl Facet for AllelesFacet {
    fn name(&self) -> &str { "Alleles" }
    fn as_any(&self) -> &dyn Any { self }
}
