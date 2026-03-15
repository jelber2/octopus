use std::sync::atomic::{AtomicBool, Ordering};

pub static DEBUG_MODE: AtomicBool = AtomicBool::new(false);
pub static TRACE_MODE: AtomicBool = AtomicBool::new(false);

pub fn is_debug_mode() -> bool {
    DEBUG_MODE.load(Ordering::Relaxed)
}

pub fn is_trace_mode() -> bool {
    TRACE_MODE.load(Ordering::Relaxed)
}
