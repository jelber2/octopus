// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Converted from C++ to Rust

use super::caller::{Caller, CallerEnvironment, CallerOptions, CallBuffer};
use crate::core::types::phylogeny::Phylogeny;

pub struct CellCaller {
    options: CallerOptions,
    num_cells: usize,
    phylogeny: Option<Phylogeny>,
}

impl CellCaller {
    pub fn new(options: CallerOptions, num_cells: usize) -> Self {
        CellCaller { options, num_cells, phylogeny: None }
    }

    pub fn with_phylogeny(mut self, phylogeny: Phylogeny) -> Self {
        self.phylogeny = Some(phylogeny);
        self
    }
}

impl Caller for CellCaller {
    fn name(&self) -> &str { "cell" }

    fn call_variants(&self, env: &CallerEnvironment) -> Result<CallBuffer, String> {
        Ok(CallBuffer::new())
    }
}
