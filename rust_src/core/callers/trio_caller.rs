// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Converted from C++ to Rust

use super::caller::{Caller, CallerEnvironment, CallerOptions, CallBuffer};
use crate::basics::trio::Trio;

pub struct TrioCaller {
    options: CallerOptions,
    trio: Trio,
}

impl TrioCaller {
    pub fn new(options: CallerOptions, trio: Trio) -> Self {
        TrioCaller { options, trio }
    }
}

impl Caller for TrioCaller {
    fn name(&self) -> &str { "trio" }

    fn call_variants(&self, env: &CallerEnvironment) -> Result<CallBuffer, String> {
        Ok(CallBuffer::new())
    }
}
