// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Converted from C++ to Rust

use super::caller::{Caller, CallerEnvironment, CallerOptions, CallBuffer};

pub struct CancerCaller {
    options: CallerOptions,
    normal_sample: Option<String>,
    somatic_mutation_rate: f64,
}

impl CancerCaller {
    pub fn new(options: CallerOptions, normal_sample: Option<String>, somatic_mutation_rate: f64) -> Self {
        CancerCaller { options, normal_sample, somatic_mutation_rate }
    }
}

impl Caller for CancerCaller {
    fn name(&self) -> &str { "cancer" }

    fn call_variants(&self, env: &CallerEnvironment) -> Result<CallBuffer, String> {
        Ok(CallBuffer::new())
    }
}
