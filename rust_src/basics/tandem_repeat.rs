// Converted from C++ to Rust

use super::genomic_region::GenomicRegion;

#[derive(Debug, Clone)]
pub struct TandemRepeat {
    region: GenomicRegion,
    period: usize,
    sequence: String,
}

impl TandemRepeat {
    pub fn new(region: GenomicRegion, period: usize, sequence: String) -> Self {
        TandemRepeat { region, period, sequence }
    }

    pub fn region(&self) -> &GenomicRegion { &self.region }
    pub fn period(&self) -> usize { self.period }
    pub fn sequence(&self) -> &str { &self.sequence }
}
