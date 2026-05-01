use std::{
    ffi::{CStr, CString},
    os::raw::c_char,
};

use crate::{RunFileRequest, RunFileResult, run_file_request, types::FfiResult};

#[unsafe(no_mangle)]
/// Executes a bioscript run request encoded as a UTF-8 JSON C string.
///
/// # Safety
///
/// `request_json` must either be null or point to a valid, NUL-terminated C
/// string that remains alive for the duration of this call.
pub unsafe extern "C" fn bioscript_run_file_json(request_json: *const c_char) -> *mut c_char {
    let response = unsafe {
        if request_json.is_null() {
            FfiResult::<RunFileResult> {
                ok: false,
                value: None,
                error: Some("request_json was null".to_owned()),
            }
        } else {
            parse_and_run_request(request_json)
        }
    };

    encode_response(&response)
}

unsafe fn parse_and_run_request(request_json: *const c_char) -> FfiResult<RunFileResult> {
    match unsafe { CStr::from_ptr(request_json) }.to_str() {
        Ok(value) => match serde_json::from_str::<RunFileRequest>(value) {
            Ok(request) => match run_file_request(request) {
                Ok(result) => FfiResult {
                    ok: true,
                    value: Some(result),
                    error: None,
                },
                Err(error) => FfiResult::<RunFileResult> {
                    ok: false,
                    value: None,
                    error: Some(error),
                },
            },
            Err(error) => FfiResult::<RunFileResult> {
                ok: false,
                value: None,
                error: Some(format!("invalid request JSON: {error}")),
            },
        },
        Err(error) => FfiResult::<RunFileResult> {
            ok: false,
            value: None,
            error: Some(format!("request_json was not valid UTF-8: {error}")),
        },
    }
}

fn encode_response(response: &FfiResult<RunFileResult>) -> *mut c_char {
    match serde_json::to_string(response) {
        Ok(json) => match CString::new(json) {
            Ok(value) => value.into_raw(),
            Err(_) => std::ptr::null_mut(),
        },
        Err(_) => std::ptr::null_mut(),
    }
}

#[unsafe(no_mangle)]
/// Frees a string previously returned by [`bioscript_run_file_json`].
///
/// # Safety
///
/// `ptr` must be null or a pointer returned by [`CString::into_raw`] from this
/// library, and it must not be freed more than once.
pub unsafe extern "C" fn bioscript_free_string(ptr: *mut c_char) {
    if !ptr.is_null() {
        unsafe {
            let _ = CString::from_raw(ptr);
        }
    }
}
