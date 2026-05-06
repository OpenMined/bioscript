use std::{
    collections::BTreeSet,
    fs,
    path::{Path, PathBuf},
};

const MAX_PRODUCTION_LINES: usize = 500;

#[test]
fn production_rust_files_stay_under_size_limit() {
    let repo_root = repo_root();
    let documented_backlog = documented_refactor_backlog(&repo_root);
    let mut actual_oversized = BTreeSet::new();
    let mut failures = Vec::new();

    for file in production_rust_files(&repo_root) {
        let relative = relative_slash_path(&repo_root, &file);
        let source = fs::read_to_string(&file)
            .unwrap_or_else(|err| panic!("failed to read {relative}: {err}"));
        let line_count = production_line_count(&source);

        if line_count > MAX_PRODUCTION_LINES {
            actual_oversized.insert(relative.clone());

            if !documented_backlog.contains(&relative) {
                failures.push(format!(
                    "{relative} has {line_count} production lines; split it or add it to AGENTS.md"
                ));
            }
        }
    }

    for documented in documented_backlog.difference(&actual_oversized) {
        failures.push(format!(
            "{documented} is listed in AGENTS.md but is no longer above {MAX_PRODUCTION_LINES} production lines"
        ));
    }

    assert!(
        failures.is_empty(),
        "production source size guard failed:\n{}",
        failures.join("\n")
    );
}

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(2)
        .expect("bioscript-core should live under rust/")
        .to_path_buf()
}

fn production_rust_files(repo_root: &Path) -> Vec<PathBuf> {
    let rust_dir = repo_root.join("rust");
    let mut files = Vec::new();

    for package in fs::read_dir(&rust_dir).expect("failed to read rust workspace directory") {
        let package = package.expect("failed to read rust workspace entry");
        let package_name = package.file_name();
        let package_name = package_name.to_string_lossy();

        if !package_name.starts_with("bioscript-") {
            continue;
        }

        let src_dir = package.path().join("src");
        if src_dir.is_dir() {
            collect_rust_files(&src_dir, &mut files);
        }
    }

    files.sort();
    files
}

fn documented_refactor_backlog(repo_root: &Path) -> BTreeSet<String> {
    let agents_path = repo_root.join("AGENTS.md");
    let agents = fs::read_to_string(&agents_path)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", agents_path.display()));
    let mut in_backlog = false;
    let mut paths = BTreeSet::new();

    for line in agents.lines() {
        if line.starts_with("## ") {
            in_backlog = line == "## Current Refactor Backlog";
            continue;
        }

        if !in_backlog {
            continue;
        }

        if let Some(path) = markdown_code_span(line) {
            paths.insert(path.to_owned());
        }
    }

    paths
}

fn markdown_code_span(line: &str) -> Option<&str> {
    let start = line.find('`')? + 1;
    let end = line[start..].find('`')? + start;
    Some(&line[start..end])
}

fn collect_rust_files(dir: &Path, files: &mut Vec<PathBuf>) {
    for entry in
        fs::read_dir(dir).unwrap_or_else(|err| panic!("failed to read {}: {err}", dir.display()))
    {
        let entry =
            entry.unwrap_or_else(|err| panic!("failed to read entry in {}: {err}", dir.display()));
        let path = entry.path();

        if path.is_dir() {
            collect_rust_files(&path, files);
        } else if path.extension().is_some_and(|ext| ext == "rs") {
            files.push(path);
        }
    }
}

fn relative_slash_path(repo_root: &Path, file: &Path) -> String {
    file.strip_prefix(repo_root)
        .expect("file should be inside repository")
        .components()
        .map(|component| component.as_os_str().to_string_lossy())
        .collect::<Vec<_>>()
        .join("/")
}

fn production_line_count(source: &str) -> usize {
    let mut count = 0;
    let mut pending_cfg_test = false;
    let mut skipped_brace_depth = None;
    let mut brace_depth = 0usize;

    for line in source.lines() {
        let trimmed = line.trim_start();

        if let Some(target_depth) = skipped_brace_depth {
            update_brace_depth(line, &mut brace_depth);
            if brace_depth < target_depth {
                skipped_brace_depth = None;
            }
            continue;
        }

        if trimmed.starts_with("#[cfg(test)]") {
            pending_cfg_test = true;
            continue;
        }

        let starts_test_module = pending_cfg_test && trimmed.starts_with("mod ");
        pending_cfg_test = pending_cfg_test && trimmed.starts_with("#[");

        if starts_test_module {
            update_brace_depth(line, &mut brace_depth);
            skipped_brace_depth = Some(brace_depth);
            continue;
        }

        count += 1;
        update_brace_depth(line, &mut brace_depth);
    }

    count
}

fn update_brace_depth(line: &str, brace_depth: &mut usize) {
    for byte in line.bytes() {
        match byte {
            b'{' => *brace_depth += 1,
            b'}' => *brace_depth = brace_depth.saturating_sub(1),
            _ => {}
        }
    }
}
