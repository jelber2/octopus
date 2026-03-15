// Converted from C++ to Rust

use crate::basics::aligned_read::AlignedRead;

pub trait ReadTransformer: Send + Sync {
    fn transform(&self, read: &mut AlignedRead);
    fn name(&self) -> &str;
}
