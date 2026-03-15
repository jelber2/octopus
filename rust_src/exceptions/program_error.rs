// Converted from C++ to Rust

use super::error::OctopusError;

pub fn program_error(msg: impl Into<String>) -> OctopusError {
    OctopusError::ProgramError(msg.into())
}

pub fn unimplemented_feature_error(msg: impl Into<String>) -> OctopusError {
    OctopusError::UnimplementedFeatureError(msg.into())
}
