mod c_api;
mod limits;
mod run_file;
mod types;
mod variant_yaml;

pub use c_api::{bioscript_free_string, bioscript_run_file_json};
pub use run_file::run_file_request;
pub use types::{
    RunFileRequest, RunFileResult, RunVariantYamlRequest, RunVariantYamlResult,
    VariantObservationResult,
};
pub use variant_yaml::run_variant_yaml_request;

#[cfg(target_os = "android")]
pub mod android {
    use crate::{RunFileRequest, run_file_request};
    use jni::JNIEnv;
    use jni::objects::{JClass, JString};

    #[unsafe(no_mangle)]
    pub unsafe extern "system" fn Java_expo_modules_bioscript_ExpoBioscriptNativeBridge_runFileNative<
        'local,
    >(
        mut env: JNIEnv<'local>,
        _class: JClass<'local>,
        request_json: JString<'local>,
    ) -> JString<'local> {
        let request_string: String = match env.get_string(&request_json) {
            Ok(value) => value.into(),
            Err(error) => {
                return env
                    .new_string(
                        serde_json::json!({
                            "ok": false,
                            "error": format!("failed to read request json from JVM: {error}")
                        })
                        .to_string(),
                    )
                    .expect("jni new_string should succeed");
            }
        };

        let response = match serde_json::from_str::<RunFileRequest>(&request_string) {
            Ok(request) => match run_file_request(request) {
                Ok(value) => serde_json::json!({ "ok": true, "value": value }).to_string(),
                Err(error) => serde_json::json!({ "ok": false, "error": error }).to_string(),
            },
            Err(error) => {
                serde_json::json!({ "ok": false, "error": format!("invalid request JSON: {error}") })
                    .to_string()
            }
        };

        env.new_string(response)
            .expect("jni new_string should succeed")
    }
}
