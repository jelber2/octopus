// Converted from C++ to Rust

use std::fmt;

#[derive(Debug)]
pub enum OctopusError {
    UserError(String),
    ProgramError(String),
    SystemError(String),
    FileOpenError { path: String, reason: String },
    MissingFileError(String),
    MalformedFileError { path: String, reason: String },
    UnwritableFileError(String),
    MissingIndexError(String),
    UnimplementedFeatureError(String),
}

impl fmt::Display for OctopusError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OctopusError::UserError(msg) => write!(f, "User error: {}", msg),
            OctopusError::ProgramError(msg) => write!(f, "Program error: {}", msg),
            OctopusError::SystemError(msg) => write!(f, "System error: {}", msg),
            OctopusError::FileOpenError { path, reason } => write!(f, "Could not open file {}: {}", path, reason),
            OctopusError::MissingFileError(path) => write!(f, "Missing file: {}", path),
            OctopusError::MalformedFileError { path, reason } => write!(f, "Malformed file {}: {}", path, reason),
            OctopusError::UnwritableFileError(path) => write!(f, "Cannot write file: {}", path),
            OctopusError::MissingIndexError(path) => write!(f, "Missing index for: {}", path),
            OctopusError::UnimplementedFeatureError(msg) => write!(f, "Unimplemented feature: {}", msg),
        }
    }
}

impl std::error::Error for OctopusError {}
