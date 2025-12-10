"""
BioVault project integration for BioScript.

This module provides an API for creating, reading, and managing BioVault projects
from BioScript notebooks and scripts.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Union

import yaml


class TemplateType(str, Enum):
    """Available BioVault templates."""

    DYNAMIC_NEXTFLOW = "dynamic-nextflow"


class TypeExpr(str, Enum):
    """Type expressions supported by BioVault."""

    STRING = "String"
    BOOL = "Bool"
    FILE = "File"
    DIRECTORY = "Directory"
    PARTICIPANT_SHEET = "ParticipantSheet"
    GENOTYPE_RECORD = "GenotypeRecord"
    BIOVAULT_CONTEXT = "BiovaultContext"
    DATASET = "Dataset"  # Remote dataset via beaver Twin
    DATASET_TWIN = "DatasetTwin"  # Twin object from dataset

    @staticmethod
    def list_of(inner: str) -> str:
        """Create a List type expression."""
        return f"List[{inner}]"

    @staticmethod
    def map_of(value_type: str) -> str:
        """Create a Map type expression."""
        return f"Map[String, {value_type}]"

    @staticmethod
    def optional(inner: str) -> str:
        """Create an optional type expression."""
        return f"{inner}?"


@dataclass
class Parameter:
    """A project parameter definition."""

    name: str
    type: str
    description: str
    default: Optional[Any] = None
    advanced: bool = False

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for YAML serialization."""
        d = {
            "name": self.name,
            "type": self.type,
            "description": self.description,
        }
        if self.default is not None:
            d["default"] = self.default
        if self.advanced:
            d["advanced"] = self.advanced
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> Parameter:
        """Create from dictionary."""
        return cls(
            name=d["name"],
            type=d["type"],
            description=d["description"],
            default=d.get("default"),
            advanced=d.get("advanced", False),
        )


@dataclass
class Input:
    """A project input specification."""

    name: str
    type: str
    description: str
    format: Optional[str] = None
    path: Optional[str] = None
    mapping: Optional[Dict[str, str]] = None
    cli_flag: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for YAML serialization."""
        d = {
            "name": self.name,
            "type": self.type,
            "description": self.description,
        }
        if self.format:
            d["format"] = self.format
        if self.path:
            d["path"] = self.path
        if self.mapping:
            d["mapping"] = self.mapping
        if self.cli_flag:
            d["cli_flag"] = self.cli_flag
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> Input:
        """Create from dictionary."""
        return cls(
            name=d["name"],
            type=d["type"],
            description=d["description"],
            format=d.get("format"),
            path=d.get("path"),
            mapping=d.get("mapping"),
            cli_flag=d.get("cli_flag"),
        )


@dataclass
class DatasetInput:
    """A dataset input specification for loading remote datasets via beaver.

    Datasets are loaded using beaver's DatasetRegistry and returned as Twin objects
    for privacy-preserving analysis.

    Example:
        >>> dataset = DatasetInput(
        ...     name="patient_genomics",
        ...     owner="data_owner@example.com",
        ...     dataset_name="single_cell_rnaseq",
        ...     asset_key="sc_data",  # Optional: specific asset within dataset
        ...     description="Patient single-cell RNA-seq data"
        ... )
    """

    name: str
    owner: str
    dataset_name: str
    description: str
    asset_key: Optional[str] = None  # Specific asset within the dataset

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for YAML serialization."""
        d = {
            "name": self.name,
            "type": "Dataset",
            "owner": self.owner,
            "dataset_name": self.dataset_name,
            "description": self.description,
        }
        if self.asset_key:
            d["asset_key"] = self.asset_key
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> DatasetInput:
        """Create from dictionary."""
        return cls(
            name=d["name"],
            owner=d["owner"],
            dataset_name=d["dataset_name"],
            description=d["description"],
            asset_key=d.get("asset_key"),
        )

    def load(self, context: Any = None) -> Any:
        """
        Load the dataset using beaver's DatasetRegistry.

        Args:
            context: Optional BeaverContext. If not provided, uses active context.

        Returns:
            Twin object for the dataset/asset, or Dataset object if no asset_key.

        Example:
            >>> ds_input = DatasetInput(...)
            >>> twin = ds_input.load()
            >>> # Use twin.mock for public mock data
            >>> # Use twin.private for actual private data (requires authorization)
        """
        try:
            import beaver

            ctx = context or beaver.ctx()
            dataset = ctx.datasets[self.owner][self.dataset_name]

            if self.asset_key:
                return dataset[self.asset_key]
            return dataset
        except ImportError:
            raise ImportError(
                "beaver package is required for dataset inputs. "
                "Install with: pip install biovault-beaver"
            ) from None


@dataclass
class Output:
    """A project output specification."""

    name: str
    type: str
    description: str
    format: Optional[str] = None
    path: Optional[str] = None
    cli_flag: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for YAML serialization."""
        d = {
            "name": self.name,
            "type": self.type,
            "description": self.description,
        }
        if self.format:
            d["format"] = self.format
        if self.path:
            d["path"] = self.path
        if self.cli_flag:
            d["cli_flag"] = self.cli_flag
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> Output:
        """Create from dictionary."""
        return cls(
            name=d["name"],
            type=d["type"],
            description=d["description"],
            format=d.get("format"),
            path=d.get("path"),
            cli_flag=d.get("cli_flag"),
        )


@dataclass
class ProcessDefinition:
    """Definition of a workflow process."""

    name: str
    script: str
    container: Optional[str] = None
    kind: str = "bioscript"

    def to_dict(self) -> Dict[str, Any]:
        data = {
            "name": self.name,
            "script": self.script,
        }
        if self.container:
            data["container"] = self.container
        if self.kind and self.kind != "bioscript":
            data["kind"] = self.kind
        return data

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> ProcessDefinition:
        return cls(
            name=data["name"],
            script=data["script"],
            container=data.get("container"),
            kind=data.get("kind", "bioscript"),
        )


@dataclass
class SQLStore:
    """Configuration for storing results in a SQL table."""

    source: str
    table_name: str
    destination: str = "SQL()"
    participant_column: Optional[str] = None
    key_column: Optional[str] = None  # backwards compatibility

    def __post_init__(self) -> None:
        if self.participant_column is None and self.key_column is not None:
            self.participant_column = self.key_column

    def to_dict(self) -> Dict[str, Any]:
        data: Dict[str, Any] = {
            "kind": "sql",
            "destination": self.destination,
            "source": self.source,
            "table_name": self.table_name,
        }
        if self.participant_column:
            data["participant_column"] = self.participant_column
        return data

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> SQLStore:
        return cls(
            source=data["source"],
            table_name=data["table_name"],
            destination=data.get("destination", "SQL()"),
            participant_column=data.get("participant_column")
            or data.get("participant_id_column")
            or data.get("key_column"),
        )


StoreConfig = Union[SQLStore, Dict[str, Any]]


def _store_to_dict(store: StoreConfig) -> Dict[str, Any]:
    if isinstance(store, SQLStore):
        return store.to_dict()
    if isinstance(store, dict):
        return store
    raise TypeError(f"Unsupported store configuration type: {type(store)!r}")


def _store_from_dict(data: Dict[str, Any]) -> StoreConfig:
    if data.get("kind") == "sql":
        return SQLStore.from_dict(data)
    return data


@dataclass
class PipelineStep:
    """A single step within a BioVault pipeline."""

    step_id: str
    uses: str
    with_args: Optional[Mapping[str, str]] = None
    publish: Optional[Mapping[str, str]] = None
    store: Optional[Mapping[str, StoreConfig]] = None

    def to_dict(self) -> Dict[str, Any]:
        data: Dict[str, Any] = {
            "id": self.step_id,
            "uses": self.uses,
        }
        if self.with_args:
            data["with"] = dict(self.with_args)
        if self.publish:
            data["publish"] = dict(self.publish)
        if self.store:
            data["store"] = {name: _store_to_dict(config) for name, config in self.store.items()}
        return data

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> PipelineStep:
        store_data = data.get("store")
        store: Optional[Dict[str, StoreConfig]] = None
        if store_data:
            store = {name: _store_from_dict(cfg) for name, cfg in store_data.items()}
        return cls(
            step_id=data["id"],
            uses=data["uses"],
            with_args=data.get("with"),
            publish=data.get("publish"),
            store=store,
        )


def _coerce_pipeline_step(value: Union[PipelineStep, Dict[str, Any]]) -> PipelineStep:
    if isinstance(value, PipelineStep):
        return value
    if isinstance(value, dict):
        return PipelineStep.from_dict(value)
    raise TypeError(f"Pipeline step must be a dict or PipelineStep, got {type(value)!r}")


@dataclass
class BioVaultPipeline:
    """Representation of a BioVault pipeline.yaml."""

    name: str
    inputs: Dict[str, str] = field(default_factory=dict)
    steps: List[PipelineStep] = field(default_factory=list)
    version: str = "0.1.0"

    _pipeline_dir: Optional[Path] = field(default=None, init=False, repr=False)

    def to_yaml(self) -> str:
        data: Dict[str, Any] = {"name": self.name, "version": self.version}
        if self.inputs:
            data["inputs"] = self.inputs
        data["steps"] = [step.to_dict() for step in self.steps]
        return yaml.dump(
            data,
            default_flow_style=False,
            sort_keys=False,
            width=4096,
            allow_unicode=True,
        )

    def save(self, path: Union[str, Path]) -> Path:
        target_dir = Path(path)
        target_dir.mkdir(parents=True, exist_ok=True)
        yaml_path = target_dir / "pipeline.yaml"
        yaml_path.write_text(self.to_yaml())
        self._pipeline_dir = target_dir
        return yaml_path


def _get_local_bioscript_version() -> Optional[str]:
    """Return the version of the locally imported bioscript package, if available."""

    try:
        import bioscript  # type: ignore

        version = getattr(bioscript, "__version__", None)
        if not version:
            return None
        version_str = str(version).strip()
        if version_str.startswith("v"):
            version_str = version_str[1:]
        return version_str or None
    except Exception:
        return None


def _default_docker_image(version_hint: Optional[str] = None) -> str:
    """Construct the default BioScript container image tag."""

    version = (version_hint or _get_local_bioscript_version() or "latest").strip()
    if version.startswith("v"):
        version = version[1:]
    version = version or "latest"
    return f"ghcr.io/openmined/bioscript:{version}"


@dataclass
class BioVaultProject:
    """
    A BioVault project configuration.

    This class represents the structure of a BioVault project.yaml file
    and provides methods to create, load, save, and export projects.
    """

    name: str
    author: str
    workflow: str = "workflow.nf"
    description: str = ""
    template: TemplateType = TemplateType.DYNAMIC_NEXTFLOW
    version: str = "0.1.0"
    assets: List[str] = field(default_factory=list)
    parameters: List[Parameter] = field(default_factory=list)
    inputs: List[Input] = field(default_factory=list)
    outputs: List[Output] = field(default_factory=list)
    processes: List[ProcessDefinition] = field(default_factory=list)

    # Docker configuration
    docker_image: Optional[str] = None
    docker_platform: str = "linux/amd64"

    # Runtime metadata (not saved to YAML)
    _project_dir: Optional[Path] = field(default=None, init=False, repr=False)
    _asset_sources: Dict[str, Path] = field(default_factory=dict, init=False, repr=False)
    _entrypoint: Optional[str] = field(default=None, init=False, repr=False)

    def __post_init__(self) -> None:
        if not self.docker_image:
            self.docker_image = _default_docker_image()
        if not self.docker_platform:
            self.docker_platform = "linux/amd64"

    @classmethod
    def from_yaml(cls, path: Union[str, Path]) -> BioVaultProject:
        """Load a project from a project.yaml file."""
        path = Path(path)
        yaml_path = path / "project.yaml" if path.is_dir() else path

        with open(yaml_path) as f:
            data = yaml.safe_load(f)

        project = cls(
            name=data["name"],
            author=data["author"],
            workflow=data.get("workflow", "workflow.nf"),
            description=data.get("description", ""),
            template=TemplateType.DYNAMIC_NEXTFLOW,  # Only dynamic-nextflow supported
            version=data.get("version", "0.1.0"),
            assets=data.get("assets", []),
            parameters=[Parameter.from_dict(p) for p in data.get("parameters", [])],
            inputs=[Input.from_dict(i) for i in data.get("inputs", [])],
            outputs=[Output.from_dict(o) for o in data.get("outputs", [])],
            processes=[ProcessDefinition.from_dict(p) for p in data.get("processes", [])],
        )

        docker_image = data.get("docker_image")
        docker_platform = data.get("docker_platform", "linux/amd64")
        project.set_docker_image(docker_image or _default_docker_image(), docker_platform)

        if project.processes:
            project._entrypoint = project.processes[0].script
        elif project.assets and not project._entrypoint:
            project._entrypoint = project.assets[0]

        # Store the project directory for later use
        project._project_dir = yaml_path.parent
        return project

    def to_yaml(self) -> str:
        """Convert project to YAML string."""
        data = {
            "name": self.name,
            "author": self.author,
            "workflow": self.workflow,
            "template": self.template.value,
            "version": self.version,
        }

        # Docker image and platform are hardcoded in workflow generation, not exposed in YAML
        # if self.docker_image:
        #     data["docker_image"] = self.docker_image
        # if self.docker_platform:
        #     data["docker_platform"] = self.docker_platform

        if self.assets:
            data["assets"] = self.assets

        if self.description:
            data["description"] = self.description

        if self.parameters:
            data["parameters"] = [p.to_dict() for p in self.parameters]

        if self.inputs:
            data["inputs"] = [i.to_dict() for i in self.inputs]

        if self.outputs:
            data["outputs"] = [o.to_dict() for o in self.outputs]

        return yaml.dump(
            data,
            default_flow_style=False,
            sort_keys=False,
            width=4096,
            allow_unicode=True,
        )

    def save(self, path: Optional[Union[str, Path]] = None) -> Path:
        """
        Save the project to a directory.

        Args:
            path: Directory to save to. If None, uses the project's original directory.

        Returns:
            Path to the saved project.yaml file.
        """
        if path is None and self._project_dir is None:
            raise ValueError("No path specified and project has no directory set")

        project_dir = Path(path) if path else self._project_dir
        project_dir.mkdir(parents=True, exist_ok=True)

        yaml_path = project_dir / "project.yaml"
        with open(yaml_path, "w") as f:
            f.write(self.to_yaml())

        self._project_dir = project_dir
        return yaml_path

    def add_parameter(
        self,
        name: str,
        param_type: str = "String",
        description: str = "",
        default: Optional[Any] = None,
        advanced: bool = False,
    ) -> BioVaultProject:
        """Add a parameter to the project."""
        self.parameters.append(
            Parameter(
                name=name,
                type=param_type,
                description=description,
                default=default,
                advanced=advanced,
            )
        )
        return self

    def add_input(
        self,
        name: str,
        input_type: str,
        description: str,
        format: Optional[str] = None,
        path: Optional[str] = None,
        mapping: Optional[Dict[str, str]] = None,
        cli_flag: Optional[str] = None,
    ) -> BioVaultProject:
        """Add an input to the project."""
        self.inputs.append(
            Input(
                name=name,
                type=input_type,
                description=description,
                format=format,
                path=path,
                mapping=mapping,
                cli_flag=cli_flag,
            )
        )
        return self

    def add_output(
        self,
        name: str,
        output_type: str,
        description: str,
        format: Optional[str] = None,
        path: Optional[str] = None,
        cli_flag: Optional[str] = None,
    ) -> BioVaultProject:
        """Add an output to the project."""
        self.outputs.append(
            Output(
                name=name,
                type=output_type,
                description=description,
                format=format,
                path=path,
                cli_flag=cli_flag,
            )
        )
        return self

    def add_asset(
        self, asset_path: str, *, source: Optional[Union[str, Path]] = None
    ) -> BioVaultProject:
        """Add an asset file to the project."""
        if asset_path not in self.assets:
            self.assets.append(asset_path)
        if source is not None:
            self._asset_sources[asset_path] = Path(source)
        return self

    def set_entrypoint(self, asset_path: str) -> BioVaultProject:
        """Set the primary asset used as the workflow entrypoint."""
        if asset_path not in self.assets:
            raise ValueError(
                f"Entrypoint asset '{asset_path}' must be added to the project before assignment"
            )
        self._entrypoint = asset_path
        return self

    def set_docker_image(
        self, image: Optional[str], platform: str = "linux/amd64"
    ) -> BioVaultProject:
        """Set the Docker image for the project."""
        resolved_image = image or _default_docker_image()
        self.docker_image = resolved_image
        self.docker_platform = platform
        return self

    def _generate_participant_workflow_nf(self, entrypoint: Optional[str] = None) -> str:
        """Generate workflow for List[GenotypeRecord] with participant iteration and aggregation."""

        primary_process = self.processes[0]
        container_image = primary_process.container or self.docker_image or _default_docker_image()
        workflow_script_asset = entrypoint or self._entrypoint or primary_process.script

        # Determine output pattern from outputs
        individual_pattern = None
        aggregated_path = None
        classifier_name = None

        for output_spec in self.outputs:
            if output_spec.path:
                if "{participant_id}" in output_spec.path:
                    individual_pattern = output_spec.path.replace("{participant_id}", "*")
                else:
                    aggregated_path = output_spec.path
                    # Extract classifier name from aggregated path (e.g., result_HERC2.tsv -> HERC2)
                    if aggregated_path.startswith("result_") and aggregated_path.endswith(".tsv"):
                        classifier_name = aggregated_path[7:-4]  # Remove "result_" and ".tsv"

        if not classifier_name:
            classifier_name = self.name.upper().replace("-", "_").replace(" ", "_")

        if not individual_pattern:
            individual_pattern = f"result_{classifier_name}_*.tsv"
        if not aggregated_path:
            aggregated_path = f"result_{classifier_name}.tsv"

        header_comment = f"// BioVault workflow export v{self.version}\n\n"

        workflow = f'''{header_comment}nextflow.enable.dsl=2

workflow USER {{
    take:
        context
        participants  // Channel emitting GenotypeRecord maps

    main:
        def assetsDir = context.assets_dir
        if (!assetsDir) {{
            throw new IllegalStateException("Missing assets directory in context")
        }}
        def assetsDirPath = file(assetsDir)

        // Pair the assets directory with each (participant_id, genotype_file) tuple
        def participant_work_items = participants.map {{ record ->
            tuple(
                assetsDirPath,
                record.participant_id,
                file(record.genotype_file)
            )
        }}

        // Process each participant
        def per_participant_results = {primary_process.name}(
            participant_work_items
        )

        // Aggregate all results into single file
        def aggregated = aggregate_results(
            per_participant_results.collect()
        )

    emit:
        {self.outputs[0].name if self.outputs else "classification_result"} = aggregated
}}

process {primary_process.name} {{
    container '{container_image}'
    publishDir params.results_dir, mode: 'copy', overwrite: true, pattern: '{individual_pattern}'
    tag {{ participant_id }}
    errorStrategy {{ params.nextflow.error_strategy }}
    maxRetries {{ params.nextflow.max_retries }}

    input:
        tuple path(assets_dir), val(participant_id), path(genotype_file)

    output:
        path "result_{classifier_name}_${{participant_id}}.tsv"

    script:
    def genoFileName = genotype_file.getName()
    """
    GENO_FILE=\\$(printf '%q' "${{genoFileName}}")
    bioscript classify "${{assets_dir}}/{workflow_script_asset}" --file \\$GENO_FILE --participant_id "${{participant_id}}"
    """
}}

process aggregate_results {{
    container '{container_image}'
    publishDir params.results_dir, mode: 'copy', overwrite: true

    input:
        path individual_results

    output:
        path "{aggregated_path}"

    script:
    def manifestContent = individual_results.collect {{ it.toString() }}.join('\\n') + '\\n'
    """
    cat <<'EOF' > results.list\\n${{manifestContent}}EOF
    bioscript combine --list results.list --output {aggregated_path}
    """
}}
'''
        return workflow

    def generate_workflow_nf(self, entrypoint: Optional[str] = None) -> str:
        """Generate a Nextflow workflow file for this workflow."""

        if not self.processes:
            script_candidate = entrypoint or self._entrypoint
            if not script_candidate and self.assets:
                script_candidate = self.assets[0]
            if not script_candidate:
                raise ValueError(
                    "BioVaultProject requires at least one asset to generate the workflow"
                )
            container_candidate = self.docker_image or _default_docker_image()
            default_name = (self.name or "process").replace(" ", "_")
            self.processes = [
                ProcessDefinition(
                    name=default_name,
                    script=script_candidate,
                    container=container_candidate,
                )
            ]
            if not self._entrypoint:
                self._entrypoint = script_candidate

        # Check if using List[GenotypeRecord] - requires different workflow pattern
        uses_genotype_list = any(inp.type.startswith("List[GenotypeRecord") for inp in self.inputs)

        if uses_genotype_list:
            return self._generate_participant_workflow_nf(entrypoint)

        if len(self.processes) > 1:
            raise NotImplementedError("Multiple processes per workflow are not supported yet")

        primary_process = self.processes[0]
        if primary_process.kind != "bioscript":
            raise NotImplementedError(
                f"Unsupported process kind '{primary_process.kind}' for workflow generation"
            )

        if primary_process.script not in self.assets:
            raise ValueError(
                f"Primary process script '{primary_process.script}' is not registered as an asset"
            )

        entrypoint_asset = entrypoint or self._entrypoint or primary_process.script
        if entrypoint_asset not in self.assets:
            raise ValueError(f"Entrypoint asset '{entrypoint_asset}' is not registered as an asset")

        def _identifier(name: str) -> str:
            sanitized = "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in name)
            if sanitized and sanitized[0].isdigit():
                sanitized = f"_{sanitized}"
            return sanitized or "value"

        process_symbol = _identifier(primary_process.name or self.name or "process")

        workflow_inputs = [("context", "context")]
        for input_spec in self.inputs:
            workflow_inputs.append((input_spec.name, _identifier(input_spec.name)))

        take_section = "\n".join(f"        {nf_name}" for _, nf_name in workflow_inputs)
        process_args = ["script_ch"] + [nf_name for _, nf_name in workflow_inputs[1:]]

        output_channels: List[str] = []
        emit_lines: List[str] = []
        for output_spec in self.outputs:
            sanitized_name = _identifier(output_spec.name)
            channel_name = f"{sanitized_name}_ch"
            output_channels.append(channel_name)
            emit_lines.append(f"        {sanitized_name} = {channel_name}")

        if not output_channels:
            process_assignment = None
        elif len(output_channels) == 1:
            process_assignment = f"        def {output_channels[0]} = {process_symbol}(\n"
        else:
            tuple_vars = ", ".join(output_channels)
            process_assignment = f"        def ({tuple_vars}) = {process_symbol}(\n"

        input_specs = ["        path script"]
        cli_lines: List[str] = []
        for idx, input_spec in enumerate(self.inputs):
            nf_name = _identifier(input_spec.name)
            input_type = input_spec.type.lower()
            if input_type in {"file", "directory"}:
                input_specs.append(f"        path {nf_name}")
            else:
                input_specs.append(f"        val {nf_name}")

            flag = input_spec.cli_flag or (
                "--input" if idx == 0 else f"--{input_spec.name.replace('_', '-')}"
            )
            cli_lines.append(f'        {flag} "${{{nf_name}}}"')

        output_specs: List[str] = []
        for idx, output_spec in enumerate(self.outputs):
            output_path = output_spec.path or f"{output_spec.name}.txt"
            emit_target = _identifier(output_spec.name)
            output_specs.append(f"        path '{output_path}', emit: {emit_target}")

            flag = output_spec.cli_flag or (
                "--output" if idx == 0 else f"--{output_spec.name.replace('_', '-')}"
            )
            cli_lines.append(f'        {flag} "{output_path}"')

        cli_section = "\n".join(cli_lines)
        if cli_section:
            cli_section = " \\n" + cli_section

        emit_section = "\n".join(emit_lines)

        container_image = primary_process.container or self.docker_image or _default_docker_image()
        primary_process.container = container_image
        workflow_script_asset = primary_process.script

        workflow_lines: List[str] = [
            f"// BioVault workflow export v{self.version}\n\n",
            "nextflow.enable.dsl=2\n\n",
            "workflow USER {\n",
            "    take:\n",
            take_section + "\n\n",
            "    main:\n",
            "        def assetsDir = file(context.params.assets_dir)\n",
            f'        def workflowScript = file("${{assetsDir}}/{workflow_script_asset}")\n',
            "        def script_ch = Channel.value(workflowScript)\n",
        ]

        args_block = ",\n            ".join(process_args)

        if process_assignment:
            workflow_lines.append(process_assignment)
            workflow_lines.append(f"            {args_block}\n        )\n")
        elif process_args:
            workflow_lines.append(f"        {process_symbol}(\n")
            workflow_lines.append(f"            {args_block}\n        )\n")
        else:
            workflow_lines.append(f"        {process_symbol}(script_ch)\n")

        if emit_section:
            workflow_lines.append("\n    emit:\n")
            workflow_lines.append(emit_section + "\n")

        workflow_lines.append("}\n\n")
        workflow_lines.append(f"process {process_symbol} {{\n")
        workflow_lines.append(f"    container '{container_image}'\n")
        workflow_lines.append(
            "    publishDir params.results_dir, mode: 'copy', overwrite: true\n\n"
        )
        workflow_lines.append("    input:\n")
        workflow_lines.append("\n".join(input_specs) + "\n\n")

        if output_specs:
            workflow_lines.append("    output:\n")
            workflow_lines.append("\n".join(output_specs) + "\n\n")

        workflow_lines.append("    script:\n")
        workflow_lines.append('    """\n')
        workflow_lines.append("    python3 ${script}" + cli_section + "\n")
        workflow_lines.append('    """\n')
        workflow_lines.append("}\n")

        return "".join(workflow_lines)

    def export(
        self,
        target_dir: Union[str, Path],
        notebook_path: Optional[Union[str, Path]] = None,
        classifier_name: Optional[str] = None,
    ) -> Path:
        """
        Export the project to a directory with all necessary files.

        Args:
            target_dir: Directory to export to
            notebook_path: Optional Jupyter notebook to convert and include
            classifier_name: Optional name for the classifier script

        Returns:
            Path to the exported project directory.
        """
        target = Path(target_dir)
        target.mkdir(parents=True, exist_ok=True)

        if not self.processes:
            script_candidate = self._entrypoint or (self.assets[0] if self.assets else None)
            if not script_candidate:
                raise ValueError(
                    "Cannot export workflow without at least one asset to use as entrypoint"
                )
            container_candidate = self.docker_image or _default_docker_image()
            self.processes = [
                ProcessDefinition(
                    name=(self.name or "process").replace(" ", "_"),
                    script=script_candidate,
                    container=container_candidate,
                )
            ]
            self._entrypoint = script_candidate

        if not self._entrypoint and self.processes:
            self._entrypoint = self.processes[0].script

        if not self.docker_image:
            self.set_docker_image(_default_docker_image(), self.docker_platform)

        # Save project.yaml
        self.save(target)

        # Generate and save workflow.nf
        workflow_path = target / "workflow.nf"
        with open(workflow_path, "w") as f:
            f.write(self.generate_workflow_nf(entrypoint=self._entrypoint))

        # Create assets directory
        assets_dir = target / "assets"
        assets_dir.mkdir(exist_ok=True)

        copied_assets = set()

        # Copy assets provided explicitly
        for asset, src in self._asset_sources.items():
            src_path = Path(src)
            if not src_path.exists():
                raise FileNotFoundError(f"Asset source not found: {src_path}")
            dst_path = assets_dir / asset
            dst_path.parent.mkdir(parents=True, exist_ok=True)
            if src_path.is_dir():
                shutil.copytree(src_path, dst_path, dirs_exist_ok=True)
            else:
                shutil.copy2(src_path, dst_path)
            copied_assets.add(asset)

        # Export notebook if provided
        if notebook_path:
            from .testing import export_from_notebook

            notebook_path = Path(notebook_path)
            script_name = classifier_name or notebook_path.stem + ".py"
            script_path = assets_dir / script_name

            # Export notebook to Python script
            export_from_notebook(notebook_path, script_path)

            # Add to assets if not already there
            if script_name not in self.assets:
                self.assets.append(script_name)
                self.save(target)  # Update project.yaml
            copied_assets.add(script_name)

        # Copy existing assets if project has a source directory
        if self._project_dir:
            src_assets = self._project_dir / "assets"
            if src_assets.exists():
                for asset in self.assets:
                    if asset in copied_assets:
                        continue
                    src = src_assets / asset
                    if src.exists():
                        dst = assets_dir / asset
                        # Only copy if source and destination are different
                        if src.resolve() != dst.resolve():
                            if src.is_dir():
                                shutil.copytree(src, dst, dirs_exist_ok=True)
                            else:
                                shutil.copy2(src, dst)

        return target

    def scaffold(self, target_dir: Union[str, Path]) -> Path:
        """
        Create a scaffold project structure using the BioVault CLI.

        This requires the 'bv' command to be available in the PATH.

        Args:
            target_dir: Directory to create the project in

        Returns:
            Path to the created project.
        """
        target = Path(target_dir)

        # Save project.yaml to a temporary location
        import tempfile

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            f.write(self.to_yaml())
            temp_yaml = f.name

        try:
            # Run bv project init command
            result = subprocess.run(
                ["bv", "project", "init", "--spec", temp_yaml, "--target", str(target)],
                capture_output=True,
                text=True,
            )

            if result.returncode != 0:
                raise RuntimeError(f"Failed to scaffold project: {result.stderr}")

        finally:
            # Clean up temporary file
            Path(temp_yaml).unlink(missing_ok=True)

        self._project_dir = target
        return target


def create_bioscript_project(
    name: str,
    author: str,
    description: str = "",
    docker_image: Optional[str] = None,
) -> BioVaultProject:
    """
    Create a new BioVault project configured for BioScript classifiers.

    Args:
        name: Project name (should be kebab-case)
        author: Author email
        description: Project description
        docker_image: Docker image to use for execution

    Returns:
        A configured BioVaultProject instance.
    """
    project = BioVaultProject(
        name=name,
        author=author,
        template=TemplateType.DYNAMIC_NEXTFLOW,
        description=description,
    )

    # Set Docker image
    project.set_docker_image(docker_image or _default_docker_image())

    # Add standard BioScript inputs
    project.add_input(
        name="genotype_file",
        input_type="File",
        description="TSV file containing genotype data",
        format="tsv",
    )

    # Add standard BioScript outputs
    project.add_output(
        name="classification_result",
        output_type="File",
        description=description or f"{name} classification results",
        format="tsv",
        path=f"{name}_results.tsv",
    )

    # Add common parameters
    project.add_parameter(
        name="participant_id",
        param_type="String",
        description="Participant identifier",
        default="unknown",
    )

    return project


def _coerce_parameter(value: Union[Parameter, Dict[str, Any]]) -> Parameter:
    if isinstance(value, Parameter):
        return value
    return Parameter.from_dict(value)


def _coerce_input(value: Union[Input, Dict[str, Any]]) -> Input:
    if isinstance(value, Input):
        return value
    return Input.from_dict(value)


def _coerce_output(value: Union[Output, Dict[str, Any]]) -> Output:
    if isinstance(value, Output):
        return value
    return Output.from_dict(value)


def _coerce_process(value: Union[ProcessDefinition, Dict[str, Any]]) -> ProcessDefinition:
    if isinstance(value, ProcessDefinition):
        return value
    if not isinstance(value, dict):
        raise TypeError(
            f"Process definition must be a dict or ProcessDefinition, got {type(value)}"
        )
    return ProcessDefinition.from_dict(value)


def export_workflow(
    workflow_name: str,
    author: str,
    target_dir: Union[str, Path],
    *,
    processes: Sequence[Union[ProcessDefinition, Dict[str, Any]]],
    assets: Optional[Mapping[str, Union[str, Path]]] = None,
    inputs: Optional[Sequence[Union[Input, Dict[str, Any]]]] = None,
    outputs: Optional[Sequence[Union[Output, Dict[str, Any]]]] = None,
    parameters: Optional[Sequence[Union[Parameter, Dict[str, Any]]]] = None,
    version: str = "0.1.0",
    description: str = "",
    docker_image: Optional[str] = None,
) -> BioVaultProject:
    """Export a BioScript workflow with explicit process definitions."""

    if not processes:
        raise ValueError("export_workflow requires at least one process definition")

    workflow_root = Path(target_dir) / workflow_name

    project = BioVaultProject(
        name=workflow_name,
        author=author,
        template=TemplateType.DYNAMIC_NEXTFLOW,
        version=version,
        description=description,
    )
    project.set_docker_image(docker_image)

    # Reset collections configured by the dataclass defaults
    project.parameters = []
    project.inputs = []
    project.outputs = []
    project.assets = []
    project._asset_sources.clear()
    project._entrypoint = None
    project.processes = []

    if assets:
        for dest, src in assets.items():
            dest_name = str(dest)
            src_path = Path(src)
            if not src_path.is_absolute():
                src_path = src_path.resolve()
            project.add_asset(dest_name, source=src_path)

    if parameters:
        project.parameters = [_coerce_parameter(p) for p in parameters]

    if inputs:
        project.inputs = [_coerce_input(i) for i in inputs]

    if outputs:
        project.outputs = [_coerce_output(o) for o in outputs]

    project.processes = [_coerce_process(p) for p in processes]
    if not project.processes:
        raise ValueError("export_workflow requires at least one process definition")

    for proc in project.processes:
        if proc.script not in project.assets:
            raise ValueError(f"Process script '{proc.script}' is not present in registered assets")
        if not proc.container:
            proc.container = project.docker_image

    primary_process = project.processes[0]
    project.set_entrypoint(primary_process.script)

    project.export(workflow_root)
    return project


def export_bioscript_workflow(
    script_path: Union[str, Path],
    workflow_name: str,
    author: str,
    target_dir: Union[str, Path],
    *,
    assets: Optional[
        Union[Mapping[str, Union[str, Path]], Sequence[Union[str, Path]], str, Path]
    ] = None,
    inputs: Optional[Sequence[Union[Input, Dict[str, Any]]]] = None,
    outputs: Optional[Sequence[Union[Output, Dict[str, Any]]]] = None,
    parameters: Optional[Sequence[Union[Parameter, Dict[str, Any]]]] = None,
    version: str = "0.1.0",
    description: str = "",
    docker_image: Optional[str] = None,
) -> BioVaultProject:
    """Convenience helper to export a single-script BioScript workflow."""

    script_path = Path(script_path).resolve()
    base_dir = script_path.parent

    combined_assets: Dict[str, Union[str, Path]] = {script_path.name: script_path}
    if assets:
        if isinstance(assets, Mapping):
            items = assets.items()
        else:
            if isinstance(assets, (str, Path)):
                asset_iterable: Sequence[Union[str, Path]] = [assets]
            else:
                asset_iterable = list(assets)
            items = ((Path(src).name, src) for src in asset_iterable)

        for dest, src in items:
            dest_name = str(dest)
            src_path = Path(src)
            if not src_path.is_absolute():
                src_path = (base_dir / src_path).resolve()
            if dest_name == script_path.name:
                continue
            combined_assets[dest_name] = src_path

    process_name = workflow_name.replace("-", "_")
    container = docker_image if docker_image else None
    process = ProcessDefinition(name=process_name, script=script_path.name, container=container)

    return export_workflow(
        workflow_name=workflow_name,
        author=author,
        target_dir=target_dir,
        processes=[process],
        assets=combined_assets,
        inputs=inputs,
        outputs=outputs,
        parameters=parameters,
        version=version,
        description=description,
        docker_image=docker_image,
    )


def export_bioscript_pipeline(
    pipeline_name: str,
    target_dir: Union[str, Path],
    *,
    inputs: Optional[Mapping[str, str]] = None,
    steps: Sequence[Union[PipelineStep, Dict[str, Any]]],
    version: str = "0.1.0",
) -> BioVaultPipeline:
    """Export a BioVault pipeline.yaml file."""

    if not steps:
        raise ValueError("export_bioscript_pipeline requires at least one step definition")

    pipeline_steps = [_coerce_pipeline_step(step) for step in steps]
    pipeline = BioVaultPipeline(
        name=pipeline_name,
        inputs=dict(inputs or {}),
        steps=list(pipeline_steps),
        version=version,
    )
    pipeline.save(target_dir)
    return pipeline


def export_notebook_as_project(
    notebook_path: Union[str, Path],
    project_name: str,
    author: str,
    target_dir: Union[str, Path],
    description: str = "",
    docker_image: Optional[str] = None,
) -> BioVaultProject:
    """
    Export a Jupyter notebook as a packaged BioVault workflow.

    This function:
    1. Creates a BioVault project configuration
    2. Exports the notebook as a Python script
    3. Generates the workflow.nf file
    4. Creates the complete project structure

    Args:
        notebook_path: Path to the Jupyter notebook
        project_name: Name for the BioVault project
        author: Author email
        target_dir: Directory to export the project to
        description: Project description
        docker_image: Docker image to use

    Returns:
        The created BioVaultProject instance.
    """
    # Create the project
    project = create_bioscript_project(
        name=project_name,
        author=author,
        description=description,
        docker_image=docker_image,
    )

    # Extract classifier name from notebook
    notebook_path = Path(notebook_path)
    classifier_script = notebook_path.stem + "_classifier.py"
    project.add_asset(classifier_script)
    project.set_entrypoint(classifier_script)

    process_name = project_name.replace("-", "_")
    project.processes = [
        ProcessDefinition(
            name=process_name,
            script=classifier_script,
            container=project.docker_image,
        )
    ]

    # Add BioScript library files as assets
    bioscript_assets = [
        "bioscript/__init__.py",
        "bioscript/classifier.py",
        "bioscript/reader.py",
        "bioscript/types.py",
        "bioscript/testing.py",
        "bioscript/data.py",
    ]
    for asset in bioscript_assets:
        project.add_asset(asset)

    # Export the complete project
    project.export(
        target_dir=target_dir,
        notebook_path=notebook_path,
        classifier_name=classifier_script,
    )

    # Copy BioScript library files
    target = Path(target_dir)
    bioscript_dir = target / "assets" / "bioscript"
    bioscript_dir.mkdir(exist_ok=True)

    # Copy bioscript module files
    import bioscript

    bioscript_src = Path(bioscript.__file__).parent

    for file in ["__init__.py", "classifier.py", "reader.py", "types.py", "testing.py", "data.py"]:
        src_file = bioscript_src / file
        if src_file.exists():
            shutil.copy2(src_file, bioscript_dir / file)

    return project


# Convenience functions for notebook usage
def load_project(path: Union[str, Path]) -> BioVaultProject:
    """Load a BioVault project from a directory or YAML file."""
    return BioVaultProject.from_yaml(path)


def new_project(name: str, author: str) -> BioVaultProject:
    """Create a new empty BioVault project."""
    return BioVaultProject(name=name, author=author)
