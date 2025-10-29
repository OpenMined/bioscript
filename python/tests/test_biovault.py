"""Tests for BioVault project integration."""

import tempfile
from pathlib import Path

import yaml

from bioscript.biovault import (
    BioVaultProject,
    Input,
    Output,
    Parameter,
    TemplateType,
    TypeExpr,
    create_bioscript_project,
    load_project,
    new_project,
)


def test_template_type_enum():
    """Test that only dynamic-nextflow template is available."""
    assert TemplateType.DYNAMIC_NEXTFLOW == "dynamic-nextflow"
    # Ensure no other templates exist
    assert len(list(TemplateType)) == 1


def test_type_expr_helpers():
    """Test type expression helper methods."""
    assert TypeExpr.list_of("File") == "List[File]"
    assert TypeExpr.map_of("String") == "Map[String, String]"
    assert TypeExpr.optional("Directory") == "Directory?"


def test_parameter_creation():
    """Test Parameter creation and serialization."""
    param = Parameter(
        name="threshold",
        type="String",
        description="Quality threshold",
        default="30",
        advanced=True,
    )

    d = param.to_dict()
    assert d["name"] == "threshold"
    assert d["type"] == "String"
    assert d["description"] == "Quality threshold"
    assert d["default"] == "30"
    assert d["advanced"] is True

    # Test round-trip
    param2 = Parameter.from_dict(d)
    assert param2.name == param.name
    assert param2.type == param.type
    assert param2.default == param.default


def test_input_output_creation():
    """Test Input and Output creation."""
    inp = Input(
        name="genotypes",
        type="File",
        description="Genotype data",
        format="tsv",
        path="genotypes.tsv",
    )

    out = Output(
        name="results",
        type="File",
        description="Classification results",
        format="csv",
        path="results.csv",
    )

    assert inp.to_dict()["format"] == "tsv"
    assert out.to_dict()["path"] == "results.csv"


def test_biovault_project_creation():
    """Test BioVaultProject creation with defaults."""
    project = BioVaultProject(
        name="test-project",
        author="test@example.com",
    )

    assert project.name == "test-project"
    assert project.author == "test@example.com"
    assert project.workflow == "workflow.nf"
    assert project.template == TemplateType.DYNAMIC_NEXTFLOW
    assert project.version == "0.1.0"
    assert project.docker_image == "ghcr.io/openmined/bioscript:latest"


def test_project_add_methods():
    """Test adding parameters, inputs, and outputs to a project."""
    project = new_project("test", "test@example.com")

    # Add parameter
    project.add_parameter(
        name="param1",
        param_type="Bool",
        description="Test parameter",
        default=False,
    )
    assert len(project.parameters) == 1
    assert project.parameters[0].name == "param1"

    # Add input
    project.add_input(
        name="input1",
        input_type="File",
        description="Test input",
        format="tsv",
    )
    assert len(project.inputs) == 1
    assert project.inputs[0].format == "tsv"

    # Add output
    project.add_output(
        name="output1",
        output_type="Directory",
        description="Test output",
    )
    assert len(project.outputs) == 1
    assert project.outputs[0].type == "Directory"

    # Add asset
    project.add_asset("script.py")
    assert "script.py" in project.assets


def test_project_yaml_serialization():
    """Test project serialization to YAML."""
    project = create_bioscript_project(
        name="test-classifier",
        author="test@example.com",
        description="Test project",
    )

    yaml_str = project.to_yaml()
    data = yaml.safe_load(yaml_str)

    assert data["name"] == "test-classifier"
    assert data["author"] == "test@example.com"
    assert data["template"] == "dynamic-nextflow"
    assert data["version"] == "0.1.0"

    # Check default input/output
    assert len(data["inputs"]) == 1
    assert data["inputs"][0]["name"] == "genotype_file"
    assert len(data["outputs"]) == 1
    assert data["outputs"][0]["name"] == "classification_result"


def test_project_save_and_load():
    """Test saving and loading a project."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)

        # Create and save project
        project = new_project("test", "test@example.com")
        project.add_parameter("param1", "String", "Test param")

        yaml_path = project.save(tmppath / "test-project")
        assert yaml_path.exists()
        assert yaml_path.name == "project.yaml"

        # Load project
        loaded = load_project(tmppath / "test-project")
        assert loaded.name == "test"
        assert loaded.author == "test@example.com"
        assert len(loaded.parameters) == 1
        assert loaded.parameters[0].name == "param1"


def test_workflow_generation():
    """Test Nextflow workflow generation."""
    project = create_bioscript_project(
        name="test-classifier",
        author="test@example.com",
    )
    project.add_asset("classifier.py")

    workflow = project.generate_workflow_nf()

    # Check key elements
    assert "nextflow.enable.dsl=2" in workflow
    assert "workflow USER {" in workflow
    assert "process TEST_CLASSIFIER {" in workflow
    assert "container 'ghcr.io/openmined/bioscript:latest'" in workflow  # Default container
    assert "classifier.py" in workflow


def test_project_export():
    """Test project export with notebook conversion."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)

        # Create a dummy notebook file
        notebook_path = tmppath / "test.ipynb"
        notebook_content = {
            "cells": [
                {
                    "cell_type": "code",
                    "source": "print('test')",
                    "outputs": [],
                }
            ],
            "metadata": {},
            "nbformat": 4,
            "nbformat_minor": 5,
        }

        import json

        with open(notebook_path, "w") as f:
            json.dump(notebook_content, f)

        # Create project and export (without notebook for now)
        project = create_bioscript_project(
            name="test",
            author="test@example.com",
        )
        project.add_asset("test_classifier.py")

        export_path = project.export(tmppath / "exported")

        # Check exported files
        assert (export_path / "project.yaml").exists()
        assert (export_path / "workflow.nf").exists()
        assert (export_path / "assets").is_dir()


def test_load_project_from_yaml_file():
    """Test loading a project directly from a YAML file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)

        # Create project and save
        project = new_project("test", "test@example.com")
        yaml_path = project.save(tmppath / "test-project")

        # Load from YAML file directly
        loaded = load_project(yaml_path)
        assert loaded.name == "test"

        # Load from directory
        loaded2 = load_project(tmppath / "test-project")
        assert loaded2.name == "test"


def test_docker_configuration():
    """Test Docker image configuration."""
    project = new_project("test", "test@example.com")

    # Default Docker settings
    assert project.docker_image == "ghcr.io/openmined/bioscript:latest"
    assert project.docker_platform == "linux/amd64"

    # Set custom Docker image
    project.set_docker_image("python:3.11", "linux/arm64")
    assert project.docker_image == "python:3.11"
    assert project.docker_platform == "linux/arm64"

    # Add asset before generating workflow
    project.add_asset("test_classifier.py")

    # Check it appears in workflow
    workflow = project.generate_workflow_nf()
    assert "container 'python:3.11'" in workflow
