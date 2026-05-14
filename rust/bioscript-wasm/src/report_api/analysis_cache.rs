use bioscript_core::VariantObservation;

pub(super) fn analysis_cache_observations(
    manifest_observations: &[VariantObservation],
    app_observations: &[serde_json::Value],
) -> Vec<VariantObservation> {
    manifest_observations
        .iter()
        .map(|observation| {
            let mut observation = observation.clone();
            if let Some(app_observation) = matching_app_observation(&observation, app_observations)
                && let Some(genotype_display) = app_observation
                    .get("genotype_display")
                    .and_then(serde_json::Value::as_str)
                    .filter(|value| !value.is_empty() && *value != "??")
            {
                observation.genotype = Some(genotype_display.to_owned());
            }
            observation
        })
        .collect()
}

fn matching_app_observation<'a>(
    observation: &VariantObservation,
    app_observations: &'a [serde_json::Value],
) -> Option<&'a serde_json::Value> {
    let matched_rsid = observation.matched_rsid.as_deref()?;
    app_observations.iter().find(|app_observation| {
        app_observation
            .get("rsid")
            .and_then(serde_json::Value::as_str)
            == Some(matched_rsid)
    })
}
