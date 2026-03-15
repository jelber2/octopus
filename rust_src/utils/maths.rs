// Converted from C++ to Rust

pub fn log_sum_exp(logs: &[f64]) -> f64 {
    if logs.is_empty() { return f64::NEG_INFINITY; }
    let max = logs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if max.is_infinite() { return max; }
    max + logs.iter().map(|&x| (x - max).exp()).sum::<f64>().ln()
}

pub fn log_sum_exp2(a: f64, b: f64) -> f64 {
    if a >= b {
        a + (1.0 + (b - a).exp()).ln()
    } else {
        b + (1.0 + (a - b).exp()).ln()
    }
}

pub fn count_leading_zeros(x: f64) -> u32 {
    if x == 0.0 { return 0; }
    let abs_x = x.abs();
    if abs_x >= 1.0 { return 0; }
    (-abs_x.log10().floor()) as u32
}

pub fn is_in_unit_interval(x: f64) -> bool {
    x >= 0.0 && x <= 1.0
}

pub fn probability_to_phred(p: f64) -> f64 {
    if p <= 0.0 { return f64::INFINITY; }
    if p >= 1.0 { return 0.0; }
    -10.0 * p.log10()
}

pub fn phred_to_probability(phred: f64) -> f64 {
    10_f64.powf(-phred / 10.0)
}

pub fn round_to_nearest(x: f64, precision: u32) -> f64 {
    let factor = 10_f64.powi(precision as i32);
    (x * factor).round() / factor
}
