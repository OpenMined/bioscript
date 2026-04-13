# Assay Loader Contract

This document defines the shared loader contract for `bioscript:assay` packages.

The goal is to keep assay loading:

- package-driven
- validation-first
- explicit about runtime support
- reusable across clients such as:
  - `tools/test_assays.py`
  - BioVault / other apps

## Inputs

The loader consumes one assay package root directory. Today that means:

- `assay.yaml`
- `catalogue.yaml`
- one or more variant YAML files
- optionally a referenced script file for `implementation.kind: script`

## Validation Responsibilities

The shared loader is responsible for validating:

- `assay.yaml` declares `schema: bioscript:assay`
- `assay.yaml` includes `assay_id`
- `implementation.kind` is supported
- `implementation.path` exists for script assays
- `catalogue.yaml` declares `schema: bioscript:catalogue`
- catalogue entries are well-formed
- referenced variant files exist
- variant files are `bioscript:variant` records
- variant schema requirements are satisfied

Validation should fail fast and stop package loading.

## Runtime Capability Responsibilities

After structural validation, the loader must ask the runtime capability layer whether each variant is executable.

That capability truth is defined by:

- `runtime_capabilities.json`
- consumed by:
  - `runtime_capabilities.py`
  - app-side mirrors of the same spec

The loader must not invent its own support rules outside that shared capability source.

## Output Contract

The loader returns one `AssayPackage` with:

- package metadata
- implementation metadata
- runnable entries
- unsupported entries
- explicit unsupported reasons

Current Python representation:

```python
@dataclass
class VariantDefinition:
    name: str
    fields: dict[str, Any]
    source: str


@dataclass
class RunnableVariantEntry:
    variant: VariantDefinition


@dataclass
class UnsupportedVariantEntry:
    reason: str
    variant: VariantDefinition


@dataclass
class AssayPackage:
    implementation_kind: str
    manifest: dict[str, Any]
    root: Path
    runnable_variants: list[RunnableVariantEntry]
    script_path: Path | None
    unsupported_variants: list[UnsupportedVariantEntry]
```

Compatibility helper:

- `AssayPackage.variants`
  returns only runnable `VariantDefinition` objects for legacy runner code

That compatibility property should be treated as transitional.

## Semantics

### Runnable entries

Runnable entries are package members that:

- are structurally valid
- are supported by the current BioScript runtime capability spec

Only runnable entries should be executed by runtime-facing clients.

### Unsupported entries

Unsupported entries are package members that:

- are structurally valid enough to load
- are not executable under current runtime capability rules

Unsupported entries must carry a machine-stable human-readable reason, for example:

- `deletion missing deletion_length`
- `insertions not yet supported by bioscript runtime`
- `indels not yet supported by bioscript runtime`
- `unsupported variant kind: other`

Clients should surface these entries clearly instead of silently dropping them.

## Client Expectations

### CLI runner

`tools/test_assays.py` should:

- trust the loader for runnable vs unsupported decisions
- execute runnable entries only
- report unsupported entries and reasons

### App

Apps should mirror the same contract shape:

- assay metadata
- runnable members
- unsupported members
- reasons

Apps should not guess runtime support independently of the shared capability spec.

## Current Gap

Today the contract is aligned across:

- shared Python loader
- shared capability spec
- app-side capability spec consumer

But the app still parses package YAML in TypeScript rather than consuming a single implementation of the loader.

That is acceptable for now as long as:

- the output contract stays aligned
- runtime capability logic continues to come from the shared spec

## Next Evolution

The next refinement should be one of:

1. keep dual implementations, but freeze and test the contract more rigorously
2. introduce a shared compiled/package export format so clients stop re-parsing independently

Until then, this document is the source of truth for loader behavior.

## Shared Compiled Artifact

The first shared compiled artifact now exists in Python:

- `assay_compiled.py`

It exports a normalized serializable document with:

- assay metadata
- implementation metadata
- compatibility/privacy metadata
- runnable variants
- unsupported variants with reasons

Current envelope:

```json
{
  "schema": "bioscript:assay-compiled",
  "version": "1.0",
  "assay": {},
  "implementation": {},
  "ui": {},
  "compatibility": {},
  "privacy": {},
  "runnable_variants": [],
  "unsupported_variants": []
}
```

This is the preferred next contract for clients that should not re-interpret raw YAML directly.

### Generation workflow

Assay authors should not hand-maintain `assay.compiled.yaml`.

The source of truth remains:

- `assay.yaml`
- `catalogue.yaml`
- variant YAML files

The compiled artifact must be generated from those files before publishing a package for clients such as BioVault.

Generate one package:

```bash
python3 tools/generate_assay_compiled.py assays/traits/HERC2
```

Generate every package under a tree:

```bash
python3 tools/generate_assay_compiled.py assays
```

Check that committed compiled artifacts are up to date:

```bash
python3 tools/generate_assay_compiled.py --check assays
```

Published packages intended for thin clients should include:

- runtime/source files needed for execution
- generated `assay.compiled.yaml`
