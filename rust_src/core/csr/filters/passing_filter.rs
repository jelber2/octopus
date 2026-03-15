// Converted from C++ to Rust

use crate::io::variant::vcf_record::VcfRecord;
use super::variant_call_filter::{VariantCallFilter, FilterResult};
use super::super::facets::facet::FacetWrapper;

pub struct PassingFilter;

impl VariantCallFilter for PassingFilter {
    fn name(&self) -> &str { "PASS" }

    fn filter(&self, _record: &mut VcfRecord, _facets: &[FacetWrapper]) -> FilterResult {
        FilterResult::Pass
    }
}
