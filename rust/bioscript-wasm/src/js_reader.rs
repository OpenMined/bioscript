//! A `Read + Seek` shim backed by a JS `readAt(offset, length) -> Uint8Array`
//! callback. The host JS (Node or browser) owns the file handle; we ask for
//! byte ranges on demand so a 20 GB CRAM never needs to load into memory.
//!
//! On wasm32-unknown-unknown there are no real threads, so the `Send + Sync`
//! unsafe impls below are sound — `fasta::Repository` requires them for its
//! `Arc<RwLock<dyn Adapter + Send + Sync>>` cache.

use std::io::{self, Read, Seek, SeekFrom};

use js_sys::{Function, Uint8Array};
use wasm_bindgen::JsValue;

pub struct JsReader {
    read_at: Function,
    length: u64,
    position: u64,
    label: String,
}

impl JsReader {
    pub fn new(read_at: Function, length: u64, label: impl Into<String>) -> Self {
        Self {
            read_at,
            length,
            position: 0,
            label: label.into(),
        }
    }
}

impl Read for JsReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() || self.position >= self.length {
            return Ok(0);
        }
        let remaining = self.length - self.position;
        let want = u64::try_from(buf.len()).unwrap_or(u64::MAX).min(remaining);
        let result = self
            .read_at
            .call2(
                &JsValue::NULL,
                &JsValue::from_f64(self.position as f64),
                &JsValue::from_f64(want as f64),
            )
            .map_err(|err| {
                io::Error::other(format!(
                    "{} readAt({}, {}) threw: {:?}",
                    self.label, self.position, want, err
                ))
            })?;
        let array = Uint8Array::from(result);
        let got = array.byte_length() as usize;
        if got == 0 {
            return Ok(0);
        }
        if got > buf.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "{} readAt returned {} bytes but caller asked for {}",
                    self.label,
                    got,
                    buf.len()
                ),
            ));
        }
        array.copy_to(&mut buf[..got]);
        self.position += got as u64;
        Ok(got)
    }
}

impl Seek for JsReader {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        let new_pos = match pos {
            SeekFrom::Start(n) => n as i128,
            SeekFrom::End(n) => self.length as i128 + n as i128,
            SeekFrom::Current(n) => self.position as i128 + n as i128,
        };
        if new_pos < 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("{} seek before start of stream", self.label),
            ));
        }
        self.position = new_pos as u64;
        Ok(self.position)
    }
}

// SAFETY: wasm32-unknown-unknown is single-threaded. `JsValue`/`Function` are
// `!Send + !Sync` in the general case because shared access from multiple OS
// threads would race, but under single-threaded wasm that scenario can't
// happen. The `fasta::Repository` cache requires `Send + Sync` on its adapter
// to satisfy `Arc<RwLock<...>>` — this unsafe impl lets us satisfy that bound
// without the runtime ever actually crossing a thread boundary.
//
// Applied unconditionally: on native targets this crate is only ever built
// for type-checking (there's no real `js_sys::Function` available anyway), so
// the `Send + Sync` impls are equally unobservable there.
unsafe impl Send for JsReader {}
unsafe impl Sync for JsReader {}
