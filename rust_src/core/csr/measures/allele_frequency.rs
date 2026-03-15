// Converted from C++ to Rust

use super::measure::{Measure, MeasureValue};
use crate::io::variant::vcf_record::VcfRecord;
use crate::core::csr::facets::facet::FacetWrapper;

pub struct AlleleFrequency;

impl Measure for AlleleFrequency {
    fn name(&self) -> &str { "AF" }

    fn evaluate(&self, record: &VcfRecord, _facets: &[FacetWrapper]) -> MeasureValue {
        if let Some(af_vals) = record.info_value("AF") {
            let values: Vec<MeasureValue> = af_vals.iter()
                .filter_map(|v| v.parse::<f64>().ok().map(MeasureValue::Float))
                .collect();
            if !values.is_empty() {
                return if values.len() == 1 { values.into_iter().next().unwrap() } else { MeasureValue::Vec(values) };
            }
        }
        MeasureValue::None
    }
}
