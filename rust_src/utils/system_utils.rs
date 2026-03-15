// Converted from C++ to Rust

pub fn get_max_open_files() -> u64 {
    1024
}

pub fn get_num_processors() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}
