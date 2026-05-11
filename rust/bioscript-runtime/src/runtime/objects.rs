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

pub(crate) fn kestrel_module_object() -> MontyObject {
    empty_dataclass("KestrelModule", 26)
}

pub(crate) fn pysam_module_object() -> MontyObject {
    empty_dataclass("PysamModule", 20)
}

pub(crate) fn pyfaidx_module_object() -> MontyObject {
    empty_dataclass("PyfaidxModule", 21)
}

pub(crate) fn samtools_module_object() -> MontyObject {
    empty_dataclass("SamtoolsModule", 27)
}

pub(crate) fn vcf_module_object() -> MontyObject {
    empty_dataclass("VcfModule", 22)
}

pub(crate) fn pysam_alignment_file_object(
    path: &str,
    mode: &str,
    reference_filename: Option<&str>,
    index_filename: Option<&str>,
) -> MontyObject {
    let mut attrs = vec![
        (
            MontyObject::String("path".to_owned()),
            MontyObject::String(path.to_owned()),
        ),
        (
            MontyObject::String("mode".to_owned()),
            MontyObject::String(mode.to_owned()),
        ),
    ];
    attrs.push((
        MontyObject::String("reference_filename".to_owned()),
        reference_filename.map_or(MontyObject::None, |value| {
            MontyObject::String(value.to_owned())
        }),
    ));
    attrs.push((
        MontyObject::String("index_filename".to_owned()),
        index_filename.map_or(MontyObject::None, |value| {
            MontyObject::String(value.to_owned())
        }),
    ));
    MontyObject::Dataclass {
        name: "PysamAlignmentFile".to_owned(),
        type_id: 23,
        field_names: vec![
            "path".to_owned(),
            "mode".to_owned(),
            "reference_filename".to_owned(),
            "index_filename".to_owned(),
        ],
        attrs: attrs.into(),
        frozen: true,
    }
}

pub(crate) fn pyfaidx_fasta_object(path: &str) -> MontyObject {
    MontyObject::Dataclass {
        name: "PyfaidxFasta".to_owned(),
        type_id: 24,
        field_names: vec!["path".to_owned()],
        attrs: vec![(
            MontyObject::String("path".to_owned()),
            MontyObject::String(path.to_owned()),
        )]
        .into(),
        frozen: true,
    }
}

fn empty_dataclass(name: &str, type_id: u64) -> MontyObject {
    MontyObject::Dataclass {
        name: name.to_owned(),
        type_id,
        field_names: vec![],
        attrs: vec![].into(),
        frozen: true,
    }
}

pub(crate) fn pysam_aligned_segment_object(
    segment: &bioscript_libs::pysam::AlignedSegment,
) -> MontyObject {
    MontyObject::Dataclass {
        name: "PysamAlignedSegment".to_owned(),
        type_id: 25,
        field_names: vec![
            "query_name".to_owned(),
            "reference_name".to_owned(),
            "reference_start".to_owned(),
            "reference_end".to_owned(),
            "query_sequence".to_owned(),
            "mapping_quality".to_owned(),
            "cigarstring".to_owned(),
            "is_unmapped".to_owned(),
            "is_reverse".to_owned(),
        ],
        attrs: vec![
            optional_string_attr("query_name", segment.query_name.as_deref()),
            optional_string_attr("reference_name", segment.reference_name.as_deref()),
            optional_u64_attr("reference_start", segment.reference_start),
            optional_u64_attr("reference_end", segment.reference_end),
            optional_string_attr("query_sequence", segment.query_sequence.as_deref()),
            optional_u8_attr("mapping_quality", segment.mapping_quality),
            optional_string_attr("cigarstring", segment.cigarstring.as_deref()),
            (
                MontyObject::String("is_unmapped".to_owned()),
                MontyObject::Bool(segment.is_unmapped),
            ),
            (
                MontyObject::String("is_reverse".to_owned()),
                MontyObject::Bool(segment.is_reverse),
            ),
        ]
        .into(),
        frozen: true,
    }
}

fn optional_string_attr(name: &str, value: Option<&str>) -> (MontyObject, MontyObject) {
    (
        MontyObject::String(name.to_owned()),
        value.map_or(MontyObject::None, |value| {
            MontyObject::String(value.to_owned())
        }),
    )
}

fn optional_u64_attr(name: &str, value: Option<u64>) -> (MontyObject, MontyObject) {
    (
        MontyObject::String(name.to_owned()),
        value.map_or(MontyObject::None, |value| MontyObject::Int(value as i64)),
    )
}

fn optional_u8_attr(name: &str, value: Option<u8>) -> (MontyObject, MontyObject) {
    (
        MontyObject::String(name.to_owned()),
        value.map_or(MontyObject::None, |value| {
            MontyObject::Int(i64::from(value))
        }),
    )
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
