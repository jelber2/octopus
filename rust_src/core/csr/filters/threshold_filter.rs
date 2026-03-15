// Converted from C++ to Rust

use crate::io::variant::vcf_record::VcfRecord;
use super::variant_call_filter::{VariantCallFilter, FilterResult};
use super::super::facets::facet::FacetWrapper;

pub struct ThresholdFilter {
    info_key: String,
    threshold: f64,
    hard_filter_name: String,
    greater_than: bool,
}

impl ThresholdFilter {
    pub fn new(info_key: impl Into<String>, threshold: f64, hard_filter_name: impl Into<String>, greater_than: bool) -> Self {
        ThresholdFilter {
            info_key: info_key.into(),
            threshold,
            hard_filter_name: hard_filter_name.into(),
            greater_than,
        }
    }
}

impl VariantCallFilter for ThresholdFilter {
    fn name(&self) -> &str { &self.hard_filter_name }

    fn filter(&self, record: &mut VcfRecord, _facets: &[FacetWrapper]) -> FilterResult {
        if let Some(values) = record.info_value(&self.info_key) {
            if let Some(val_str) = values.first() {
                if let Ok(val) = val_str.parse::<f64>() {
                    let passes = if self.greater_than {
                        val >= self.threshold
                    } else {
                        val <= self.threshold
                    };
                    if passes {
                        return FilterResult::Pass;
                    } else {
                        return FilterResult::Fail(vec![self.hard_filter_name.clone()]);
                    }
                }
            }
        }
        FilterResult::Pass
    }
}
