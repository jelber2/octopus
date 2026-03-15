// Converted from C++ to Rust

pub struct DenovoModel {
    mutation_rate: f64,
}

impl DenovoModel {
    pub fn new(mutation_rate: f64) -> Self {
        DenovoModel { mutation_rate }
    }

    pub fn log_probability_denovo(&self, num_mutations: usize) -> f64 {
        if num_mutations == 0 {
            (1.0 - self.mutation_rate).max(1e-300).ln()
        } else {
            self.mutation_rate.max(1e-300).ln() * num_mutations as f64
        }
    }

    pub fn log_probability_inherited(&self, num_mutations: usize) -> f64 {
        if num_mutations == 0 { 0.0 } else { f64::NEG_INFINITY }
    }
}
