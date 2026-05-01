use bioscript_core::{VariantKind, VariantSpec};
use monty::MontyObject;

pub(crate) fn bioscript_object() -> MontyObject {
    MontyObject::Dataclass {
        name: "Bioscript".to_owned(),
        type_id: 1,
        field_names: vec![],
        attrs: vec![].into(),
        frozen: true,
    }
}

pub(crate) fn genotype_file_object(handle_id: u64) -> MontyObject {
    MontyObject::Dataclass {
        name: "GenotypeFile".to_owned(),
        type_id: 2,
        field_names: vec!["handle_id".to_owned()],
        attrs: vec![(
            MontyObject::String("handle_id".to_owned()),
            MontyObject::Int(handle_id as i64),
        )]
        .into(),
        frozen: true,
    }
}

pub(crate) fn variant_object(spec: &VariantSpec) -> MontyObject {
    let mut attrs = Vec::new();
    attrs.push((
        MontyObject::String("rsids".to_owned()),
        MontyObject::List(
            spec.rsids
                .iter()
                .cloned()
                .map(MontyObject::String)
                .collect(),
        ),
    ));
    if let Some(locus) = &spec.grch37 {
        attrs.push((
            MontyObject::String("grch37".to_owned()),
            MontyObject::String(format!("{}:{}-{}", locus.chrom, locus.start, locus.end)),
        ));
    }
    if let Some(locus) = &spec.grch38 {
        attrs.push((
            MontyObject::String("grch38".to_owned()),
            MontyObject::String(format!("{}:{}-{}", locus.chrom, locus.start, locus.end)),
        ));
    }
    if let Some(reference) = &spec.reference {
        attrs.push((
            MontyObject::String("reference".to_owned()),
            MontyObject::String(reference.clone()),
        ));
    }
    if let Some(alternate) = &spec.alternate {
        attrs.push((
            MontyObject::String("alternate".to_owned()),
            MontyObject::String(alternate.clone()),
        ));
    }
    if let Some(kind) = spec.kind {
        attrs.push((
            MontyObject::String("kind".to_owned()),
            MontyObject::String(variant_kind_name(kind).to_owned()),
        ));
    }
    if let Some(length) = spec.deletion_length {
        attrs.push((
            MontyObject::String("deletion_length".to_owned()),
            MontyObject::Int(length as i64),
        ));
    }
    if !spec.motifs.is_empty() {
        attrs.push((
            MontyObject::String("motifs".to_owned()),
            MontyObject::List(
                spec.motifs
                    .iter()
                    .cloned()
                    .map(MontyObject::String)
                    .collect(),
            ),
        ));
    }

    MontyObject::Dataclass {
        name: "Variant".to_owned(),
        type_id: 3,
        field_names: vec![
            "rsids".to_owned(),
            "grch37".to_owned(),
            "grch38".to_owned(),
            "reference".to_owned(),
            "alternate".to_owned(),
            "kind".to_owned(),
            "deletion_length".to_owned(),
            "motifs".to_owned(),
        ],
        attrs: attrs.into(),
        frozen: true,
    }
}

pub(crate) fn variant_plan_object(variants: &[VariantSpec]) -> MontyObject {
    MontyObject::Dataclass {
        name: "VariantPlan".to_owned(),
        type_id: 4,
        field_names: vec!["variants".to_owned()],
        attrs: vec![(
            MontyObject::String("variants".to_owned()),
            MontyObject::List(variants.iter().map(variant_object).collect()),
        )]
        .into(),
        frozen: true,
    }
}

pub(crate) fn variant_observation_object(
    observation: &bioscript_core::VariantObservation,
) -> MontyObject {
    let mut attrs = vec![
        (
            MontyObject::String("backend".to_owned()),
            MontyObject::String(observation.backend.clone()),
        ),
        (
            MontyObject::String("matched_rsid".to_owned()),
            match &observation.matched_rsid {
                Some(value) => MontyObject::String(value.clone()),
                None => MontyObject::None,
            },
        ),
        (
            MontyObject::String("assembly".to_owned()),
            match observation.assembly {
                Some(assembly) => MontyObject::String(match assembly {
                    bioscript_core::Assembly::Grch37 => "grch37".to_owned(),
                    bioscript_core::Assembly::Grch38 => "grch38".to_owned(),
                }),
                None => MontyObject::None,
            },
        ),
        (
            MontyObject::String("genotype".to_owned()),
            match &observation.genotype {
                Some(value) => MontyObject::String(value.clone()),
                None => MontyObject::None,
            },
        ),
        (
            MontyObject::String("ref_count".to_owned()),
            observation.ref_count.map_or(MontyObject::None, |value| {
                MontyObject::Int(i64::from(value))
            }),
        ),
        (
            MontyObject::String("alt_count".to_owned()),
            observation.alt_count.map_or(MontyObject::None, |value| {
                MontyObject::Int(i64::from(value))
            }),
        ),
        (
            MontyObject::String("depth".to_owned()),
            observation.depth.map_or(MontyObject::None, |value| {
                MontyObject::Int(i64::from(value))
            }),
        ),
        (
            MontyObject::String("decision".to_owned()),
            match &observation.decision {
                Some(value) => MontyObject::String(value.clone()),
                None => MontyObject::None,
            },
        ),
        (
            MontyObject::String("raw_counts".to_owned()),
            MontyObject::Dict(
                observation
                    .raw_counts
                    .iter()
                    .map(|(base, count)| {
                        (
                            MontyObject::String(base.clone()),
                            MontyObject::Int(i64::from(*count)),
                        )
                    })
                    .collect(),
            ),
        ),
        (
            MontyObject::String("evidence".to_owned()),
            MontyObject::List(
                observation
                    .evidence
                    .iter()
                    .cloned()
                    .map(MontyObject::String)
                    .collect(),
            ),
        ),
    ];

    MontyObject::Dataclass {
        name: "VariantObservation".to_owned(),
        type_id: 5,
        field_names: vec![
            "backend".to_owned(),
            "matched_rsid".to_owned(),
            "assembly".to_owned(),
            "genotype".to_owned(),
            "ref_count".to_owned(),
            "alt_count".to_owned(),
            "depth".to_owned(),
            "decision".to_owned(),
            "raw_counts".to_owned(),
            "evidence".to_owned(),
        ],
        attrs: attrs.drain(..).collect(),
        frozen: true,
    }
}

pub(crate) fn variant_kind_name(kind: VariantKind) -> &'static str {
    match kind {
        VariantKind::Snp => "snp",
        VariantKind::Insertion => "insertion",
        VariantKind::Deletion => "deletion",
        VariantKind::Indel => "indel",
        VariantKind::Other => "other",
    }
}
