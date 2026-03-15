pub fn init() {
    let _ = env_logger::try_init();
}

pub struct InfoLogger;
pub struct WarnLogger;
pub struct ErrorLogger;
pub struct DebugLogger;
pub struct TraceLogger;
