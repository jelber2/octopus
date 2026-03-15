// Converted from C++ to Rust

#[derive(Debug, Clone, Default)]
pub struct ReadPileup {
    depth: u32,
}

impl ReadPileup {
    pub fn new() -> Self {
        ReadPileup { depth: 0 }
    }

    pub fn depth(&self) -> u32 { self.depth }
    pub fn increment(&mut self) { self.depth += 1; }
    pub fn decrement(&mut self) { if self.depth > 0 { self.depth -= 1; } }
}
