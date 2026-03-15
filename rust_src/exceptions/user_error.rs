// Converted from C++ to Rust

use super::error::OctopusError;

pub fn user_error(msg: impl Into<String>) -> OctopusError {
    OctopusError::UserError(msg.into())
}
