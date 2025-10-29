import groovy.json.JsonOutput

nextflow.enable.dsl=2

workflow USER {
    take:
        context        // BiovaultContext
        participants   // ParticipantSheet channel providing genotype rows
        genotype_dir   // Directory path containing genotype files

    main:
        def assetsDir = file(context.params.assets_dir)
        def classifierScript = file("${assetsDir}/classify_apol1_exported.py")
        def classifier_ch = Channel.value(classifierScript)

        def output_name = 'apol1_participant_genotypes.tsv'
        def output_name_ch = Channel.value(output_name)

        def genotype_dir_ch = genotype_dir.map { it.toString() }

        def participants_all_v = participants.collect()
        def participants_json_v = participants_all_v.map { rows ->
            JsonOutput.toJson(rows ?: [])
        }

        apol1_sheet_ch = build_apol1_sheet(
            participants_json_v,
            genotype_dir_ch,
            classifier_ch,
            output_name_ch
        )

    emit:
        apol1_sheet = apol1_sheet_ch
}

process build_apol1_sheet {
    tag { 'apol1-sheet' }
    container 'ghcr.io/openmined/bioscript:latest'
    publishDir params.results_dir, mode: 'copy', pattern: output_name, overwrite: true

    input:
        val participants_json
        val genotype_dir
        path classifier_script
        val output_name

    output:
        path output_name

    script:
    """
    set -euo pipefail

    cat <<'JSON' > participants.json
${participants_json}
JSON

    export BV_GENOTYPE_DIR="${genotype_dir}"

    python3 - <<'PY'
import csv
import os
import json
import pathlib
import importlib.util

from bioscript import load_variants_tsv
from bioscript.types import MatchList

participants_path = pathlib.Path('participants.json')
participants = json.loads(participants_path.read_text(encoding='utf-8'))

base_dir = pathlib.Path(os.environ['BV_GENOTYPE_DIR']).resolve()
if not base_dir.exists():
    raise FileNotFoundError(f'Genotype directory not found: {base_dir}')
classifier_path = pathlib.Path("${classifier_script}").resolve()
output_path = pathlib.Path("${output_name}")
output_path.parent.mkdir(parents=True, exist_ok=True)

spec = importlib.util.spec_from_file_location('apol1_classifier', classifier_path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
config = getattr(module, '__bioscript__', None)
if not config:
    raise SystemExit('Classifier module missing __bioscript__ configuration')

variant_calls = config['variant_calls']
classifier = config['classifier']

base_fields = []
if participants:
    first_row = participants[0]
    if isinstance(first_row, dict):
        base_fields = list(first_row.keys())

if 'participant_id' not in base_fields:
    base_fields.insert(0, 'participant_id')
if 'apol1_genotype' not in base_fields:
    base_fields.append('apol1_genotype')

with output_path.open('w', encoding='utf-8', newline='') as fh:
    writer = csv.DictWriter(fh, fieldnames=base_fields, delimiter='\t')
    writer.writeheader()

    for idx, row in enumerate(participants):
        if not isinstance(row, dict):
            raise ValueError(f'Participant row is not a mapping: {row!r}')

        pid = row.get('participant_id') or row.get('participantID') or row.get('ParticipantID') or row.get('participant')
        if not pid:
            raise ValueError(f'Participant row missing participant_id: {row!r}')

        geno_ref = row.get('genotype_file') or row.get('genotype_file_path') or row.get('genotype_path')
        if not geno_ref:
            raise ValueError(f'Participant {pid} missing genotype file reference')

        geno_path = pathlib.Path(geno_ref)
        if not geno_path.is_absolute():
            geno_path = (base_dir / geno_path).resolve()
        if not geno_path.exists():
            raise FileNotFoundError(f'Genotype file not found for participant {pid}: {geno_path}')

        variants = load_variants_tsv(str(geno_path))
        matches = MatchList(variant_calls=variant_calls).match_rows(variants)
        genotype = str(classifier(matches))

        row_out = {key: row.get(key, '') for key in base_fields}
        row_out['participant_id'] = str(pid)
        row_out['apol1_genotype'] = genotype
        writer.writerow({key: '' if value is None else str(value) for key, value in row_out.items()})

PY
    """
}
