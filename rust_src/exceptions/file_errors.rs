// Converted from C++ to Rust

use super::error::OctopusError;

pub fn file_open_error(path: impl Into<String>, reason: impl Into<String>) -> OctopusError {
    OctopusError::FileOpenError { path: path.into(), reason: reason.into() }
}

pub fn missing_file_error(path: impl Into<String>) -> OctopusError {
    OctopusError::MissingFileError(path.into())
}

pub fn malformed_file_error(path: impl Into<String>, reason: impl Into<String>) -> OctopusError {
    OctopusError::MalformedFileError { path: path.into(), reason: reason.into() }
}

pub fn unwritable_file_error(path: impl Into<String>) -> OctopusError {
    OctopusError::UnwritableFileError(path.into())
}

pub fn missing_index_error(path: impl Into<String>) -> OctopusError {
    OctopusError::MissingIndexError(path.into())
}
