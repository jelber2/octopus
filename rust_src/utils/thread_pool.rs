// Converted from C++ to Rust

use std::sync::{Arc, Mutex};
use std::thread;

pub struct ThreadPool {
    workers: Vec<thread::JoinHandle<()>>,
    num_threads: usize,
}

impl ThreadPool {
    pub fn new(num_threads: usize) -> Self {
        ThreadPool {
            workers: Vec::new(),
            num_threads,
        }
    }

    pub fn num_threads(&self) -> usize {
        self.num_threads
    }
}
