// Converted from C++ to Rust

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct MemoryFootprint {
    bytes: u64,
}

impl MemoryFootprint {
    pub fn from_bytes(bytes: u64) -> Self { MemoryFootprint { bytes } }
    pub fn from_kb(kb: u64) -> Self { MemoryFootprint { bytes: kb * 1024 } }
    pub fn from_mb(mb: u64) -> Self { MemoryFootprint { bytes: mb * 1024 * 1024 } }
    pub fn from_gb(gb: u64) -> Self { MemoryFootprint { bytes: gb * 1024 * 1024 * 1024 } }

    pub fn bytes(&self) -> u64 { self.bytes }
    pub fn kilobytes(&self) -> f64 { self.bytes as f64 / 1024.0 }
    pub fn megabytes(&self) -> f64 { self.bytes as f64 / (1024.0 * 1024.0) }
    pub fn gigabytes(&self) -> f64 { self.bytes as f64 / (1024.0 * 1024.0 * 1024.0) }
}

impl std::fmt::Display for MemoryFootprint {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.bytes >= 1024 * 1024 * 1024 {
            write!(f, "{:.2}GB", self.gigabytes())
        } else if self.bytes >= 1024 * 1024 {
            write!(f, "{:.2}MB", self.megabytes())
        } else if self.bytes >= 1024 {
            write!(f, "{:.2}KB", self.kilobytes())
        } else {
            write!(f, "{}B", self.bytes)
        }
    }
}
