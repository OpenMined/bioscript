use std::time::Duration;

#[cfg(not(target_arch = "wasm32"))]
pub(crate) struct RuntimeInstant(std::time::Instant);

#[cfg(not(target_arch = "wasm32"))]
impl RuntimeInstant {
    pub(crate) fn now() -> Self {
        Self(std::time::Instant::now())
    }

    pub(crate) fn elapsed(&self) -> Duration {
        self.0.elapsed()
    }
}

#[cfg(target_arch = "wasm32")]
pub(crate) struct RuntimeInstant;

#[cfg(target_arch = "wasm32")]
impl RuntimeInstant {
    pub(crate) fn now() -> Self {
        Self
    }

    pub(crate) fn elapsed(&self) -> Duration {
        Duration::ZERO
    }
}
