"""Build the vendored Kestrel JAR when Apache Ant is unavailable.

Upstream Kestrel uses Ant and targets Java 7. Modern JDKs reject `-source 7`,
so this helper compiles the vendored sources with Java 8 compatibility and
packages a local `kestrel.jar` for VNtyper integration tests.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
KESTREL_ROOT = ROOT / "ports" / "vntyper" / "kestrel"
DEFAULT_OUTPUT = ROOT / "ports" / "vntyper" / "test-data" / "tools" / "kestrel" / "kestrel.jar"
MAIN_CLASS = "edu.gatech.kestrel.clui.Main"
DEPENDENCY_JARS = [
    "kanalyze.jar",
    "slf4j-api-1.7.12.jar",
    "logback-core-1.1.3.jar",
    "logback-classic-1.1.3.jar",
    "java-getopt-1.0.14.jar",
    "commons-lang3-3.4.jar",
    "xstream-1.4.5.jar",
]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    missing = missing_prerequisites()
    if missing:
        raise SystemExit("Missing prerequisites: " + ", ".join(missing))

    output = Path(args.output)
    if args.dry_run:
        print("javac " + " ".join(javac_command(Path("BUILD_CLASSES"))))
        print("jar cfm " + str(output) + " MANIFEST.MF -C BUILD_CLASSES .")
        return 0
    build_jar(output)
    return 0


def build_jar(output: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="bioscript-kestrel-build-") as temp:
        build_dir = Path(temp)
        classes_dir = build_dir / "classes"
        manifest = build_dir / "MANIFEST.MF"
        classes_dir.mkdir(parents=True)
        manifest.write_text(manifest_content(output), encoding="utf-8")
        subprocess.run(javac_command(classes_dir), check=True)
        output.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["jar", "cfm", str(output), str(manifest), "-C", str(classes_dir), "."], check=True)


def javac_command(classes_dir: Path) -> list[str]:
    return [
        "javac",
        "-source",
        "8",
        "-target",
        "8",
        "-cp",
        classpath(),
        "-d",
        str(classes_dir),
        *source_files(),
    ]


def source_files() -> list[str]:
    src_root = KESTREL_ROOT / "src"
    return [
        str(path)
        for path in sorted(src_root.rglob("*.java"))
        if "/test/" not in path.as_posix()
    ]


def classpath() -> str:
    jars = [str(KESTREL_ROOT / "lib" / name) for name in DEPENDENCY_JARS]
    return os.pathsep.join(jars)


def manifest_content(output: Path) -> str:
    return "\n".join(
        [
            "Manifest-Version: 1.0",
            manifest_attribute("Main-Class", MAIN_CLASS),
            manifest_attribute("Class-Path", manifest_classpath(output)),
            "",
        ]
    )


def manifest_attribute(name: str, value: str) -> str:
    line = f"{name}: {value}"
    if len(line) <= 70:
        return line
    lines = [line[:70]]
    rest = line[70:]
    while rest:
        lines.append(" " + rest[:69])
        rest = rest[69:]
    return "\n".join(lines)


def manifest_classpath(output: Path) -> str:
    output_parent = output.parent.resolve()
    try:
        relative_lib = (KESTREL_ROOT / "lib").resolve().relative_to(output_parent)
    except ValueError:
        relative_lib = Path(os.path.relpath((KESTREL_ROOT / "lib").resolve(), output_parent))
    return " ".join(str(relative_lib / name) for name in DEPENDENCY_JARS)


def missing_prerequisites() -> list[str]:
    missing = []
    for tool in ["javac", "jar"]:
        if shutil.which(tool) is None:
            missing.append(tool)
    for jar in DEPENDENCY_JARS:
        path = KESTREL_ROOT / "lib" / jar
        if not path.exists():
            missing.append(str(path))
    if not source_files():
        missing.append(str(KESTREL_ROOT / "src"))
    return missing


if __name__ == "__main__":
    raise SystemExit(main())
