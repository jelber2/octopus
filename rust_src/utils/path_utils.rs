// Converted from C++ to Rust

use std::path::{Path, PathBuf};

pub fn resolve_path(path: &Path, base: &Path) -> PathBuf {
    if path.is_absolute() {
        path.to_path_buf()
    } else {
        base.join(path)
    }
}

pub fn extension(path: &Path) -> Option<&str> {
    path.extension().and_then(|e| e.to_str())
}

pub fn stem(path: &Path) -> Option<&str> {
    path.file_stem().and_then(|s| s.to_str())
}

pub fn parent(path: &Path) -> Option<&Path> {
    path.parent()
}

pub fn filename(path: &Path) -> Option<&str> {
    path.file_name().and_then(|n| n.to_str())
}

pub fn append_extension(path: &Path, ext: &str) -> PathBuf {
    let mut result = path.to_path_buf();
    let old_ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
    if old_ext.is_empty() {
        result.set_extension(ext);
    } else {
        result.set_extension(format!("{}.{}", old_ext, ext));
    }
    result
}
