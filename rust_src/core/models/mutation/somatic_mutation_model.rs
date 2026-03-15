// Converted from C++ to Rust

pub struct SomaticMutationModel {
    mutation_rate: f64,
}

impl SomaticMutationModel {
    pub fn new(mutation_rate: f64) -> Self {
        SomaticMutationModel { mutation_rate }
    }

    pub fn log_probability_somatic(&self, num_mutations: usize) -> f64 {
        if num_mutations == 0 {
            (1.0 - self.mutation_rate).max(1e-300).ln()
        } else {
            self.mutation_rate.max(1e-300).ln() * num_mutations as f64
        }
    }
}
