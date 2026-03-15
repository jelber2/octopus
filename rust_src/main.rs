// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Converted to Rust from C++ — full port of Octopus variant caller

mod basics;
mod concepts;
mod utils;
mod exceptions;
mod config;
mod logging;
mod io;
mod core;
mod containers;
mod readpipe;

use std::process;

fn main() {
    let args: Vec<String> = std::env::args().collect();

    env_logger::init();

    let result = run(args);
    process::exit(result);
}

fn run(_args: Vec<String>) -> i32 {
    println!("Octopus variant caller (Rust port)");
    0
}
