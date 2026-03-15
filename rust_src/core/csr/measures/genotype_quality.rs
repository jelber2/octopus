// Converted from C++ to Rust

use super::measure::{Measure, MeasureValue};
use crate::io::variant::vcf_record::VcfRecord;
use crate::core::csr::facets::facet::FacetWrapper;

pub struct GenotypeQuality;

impl Measure for GenotypeQuality {
    fn name(&self) -> &str { "GQ" }

    fn evaluate(&self, record: &VcfRecord, _facets: &[FacetWrapper]) -> MeasureValue {
        let mut gq_values = Vec::new();
        for sample in record.samples() {
            if let Some(gq) = record.get_sample_value(sample, "GQ") {
                if let Some(v) = gq.first().and_then(|s| s.parse::<i64>().ok()) {
                    gq_values.push(MeasureValue::Int(v));
                }
            }
        }
        if gq_values.len() == 1 {
            gq_values.into_iter().next().unwrap()
        } else if !gq_values.is_empty() {
            MeasureValue::Vec(gq_values)
        } else {
            MeasureValue::None
        }
    }
}
