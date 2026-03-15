// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Converted from C++ to Rust

use super::caller::{Caller, CallerEnvironment, CallerOptions, CallBuffer};

pub struct PopulationCaller {
    options: CallerOptions,
}

impl PopulationCaller {
    pub fn new(options: CallerOptions) -> Self {
        PopulationCaller { options }
    }
}

impl Caller for PopulationCaller {
    fn name(&self) -> &str { "population" }

    fn call_variants(&self, env: &CallerEnvironment) -> Result<CallBuffer, String> {
        Ok(CallBuffer::new())
    }
}
