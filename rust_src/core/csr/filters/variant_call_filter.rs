// Converted from C++ to Rust

use crate::io::variant::vcf_record::VcfRecord;
use super::super::facets::facet::FacetWrapper;

pub trait VariantCallFilter: Send + Sync {
    fn filter(&self, record: &mut VcfRecord, facets: &[FacetWrapper]) -> FilterResult;
    fn name(&self) -> &str;
}

#[derive(Debug, Clone, PartialEq)]
pub enum FilterResult {
    Pass,
    Fail(Vec<String>),
}

impl FilterResult {
    pub fn is_pass(&self) -> bool { matches!(self, FilterResult::Pass) }
}
