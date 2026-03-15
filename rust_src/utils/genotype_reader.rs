// Converted from C++ to Rust

#[derive(Debug, Clone)]
pub struct Genotype {
    alleles: Vec<String>,
    phased: bool,
}

impl Genotype {
    pub fn new(alleles: Vec<String>, phased: bool) -> Self {
        Genotype { alleles, phased }
    }

    pub fn alleles(&self) -> &[String] { &self.alleles }
    pub fn is_phased(&self) -> bool { self.phased }
    pub fn ploidy(&self) -> usize { self.alleles.len() }
}

pub fn parse_genotype(gt: &str) -> Option<Genotype> {
    if gt.contains('|') {
        let alleles = gt.split('|').map(|a| a.to_string()).collect();
        Some(Genotype::new(alleles, true))
    } else if gt.contains('/') {
        let alleles = gt.split('/').map(|a| a.to_string()).collect();
        Some(Genotype::new(alleles, false))
    } else {
        None
    }
}
