// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Converted from C++ to Rust

use super::caller::{Caller, CallerEnvironment, CallerOptions, CallBuffer};

pub struct PolycloneCaller {
    options: CallerOptions,
    num_clones: usize,
}

impl PolycloneCaller {
    pub fn new(options: CallerOptions, num_clones: usize) -> Self {
        PolycloneCaller { options, num_clones }
    }
}

impl Caller for PolycloneCaller {
    fn name(&self) -> &str { "polyclone" }

    fn call_variants(&self, env: &CallerEnvironment) -> Result<CallBuffer, String> {
        Ok(CallBuffer::new())
    }
}
