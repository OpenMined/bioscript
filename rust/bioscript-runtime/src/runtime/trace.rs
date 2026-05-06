use bioscript_core::RuntimeError;
use monty::MontyObject;

use super::{BioscriptRuntime, args::reject_kwargs};

pub(crate) fn trace_lookup_metadata(source: &str) -> (Option<String>, Option<String>) {
    if let Some(rsid) = extract_rsid(source) {
        let url = format!("https://www.ncbi.nlm.nih.gov/snp/{rsid}");
        return (Some(rsid), Some(url));
    }

    if let Some(coord) = extract_coordinate(source) {
        let lower = source.to_ascii_lowercase();
        let host = if lower.contains("grch37") || lower.contains("hg19") {
            "https://grch37.ensembl.org"
        } else {
            "https://www.ensembl.org"
        };
        let url = format!("{host}/Homo_sapiens/Location/View?r={coord}");
        return (Some(coord), Some(url));
    }

    (None, None)
}

pub(crate) fn statement_context(lines: &[&str], line_no: usize) -> String {
    let Some(start_idx) = line_no.checked_sub(1) else {
        return String::new();
    };
    let Some(first_line) = lines.get(start_idx) else {
        return String::new();
    };

    let mut out = String::from(first_line.trim());
    let mut depth = update_nesting_depth(0, first_line);
    let mut current = start_idx + 1;

    while depth > 0 {
        let Some(line) = lines.get(current) else {
            break;
        };
        if !out.is_empty() {
            out.push(' ');
        }
        out.push_str(line.trim());
        depth = update_nesting_depth(depth, line);
        current += 1;
    }

    out
}

pub(crate) fn extract_rsid(source: &str) -> Option<String> {
    let chars: Vec<char> = source.chars().collect();
    let len = chars.len();
    let mut idx = 0;
    while idx + 2 <= len {
        if chars[idx] == 'r'
            && chars.get(idx + 1) == Some(&'s')
            && (idx == 0 || !chars[idx - 1].is_ascii_alphanumeric())
        {
            let mut end = idx + 2;
            while end < len && chars[end].is_ascii_digit() {
                end += 1;
            }
            if end > idx + 2 {
                return Some(chars[idx..end].iter().collect());
            }
        }
        idx += 1;
    }
    None
}

pub(crate) fn extract_coordinate(source: &str) -> Option<String> {
    for token in source.split(|ch: char| {
        ch.is_whitespace() || matches!(ch, '"' | '\'' | ',' | ')' | '(' | '[' | ']' | '{' | '}')
    }) {
        let cleaned = token.trim_matches(|ch: char| matches!(ch, ';'));
        let normalized = cleaned.strip_prefix("chr").unwrap_or(cleaned);
        if let Some((chrom, rest)) = normalized.split_once(':')
            && !chrom.is_empty()
            && chrom.chars().all(|ch| ch.is_ascii_alphanumeric())
        {
            if let Some((start, end)) = rest.split_once('-') {
                if start.chars().all(|ch| ch.is_ascii_digit())
                    && end.chars().all(|ch| ch.is_ascii_digit())
                {
                    return Some(format!("{chrom}:{start}-{end}"));
                }
            } else if rest.chars().all(|ch| ch.is_ascii_digit()) {
                return Some(format!("{chrom}:{rest}-{rest}"));
            }
        }
    }
    None
}

pub(crate) fn host_trace(
    runtime: &BioscriptRuntime,
    args: &[MontyObject],
    kwargs: &[(MontyObject, MontyObject)],
) -> Result<MontyObject, RuntimeError> {
    reject_kwargs(kwargs, "__bioscript_trace__")?;
    if let Some(MontyObject::Int(v)) = args.first() {
        runtime
            .state
            .trace_lines
            .lock()
            .expect("trace mutex poisoned")
            .push(*v as usize);
    }
    Ok(MontyObject::None)
}

pub(crate) fn instrument_source(code: &str) -> String {
    let mut out = Vec::new();
    let mut nesting_depth = 0usize;
    let mut pending_backslash = false;
    for (idx, line) in code.lines().enumerate() {
        let line_no = idx + 1;
        let trimmed = line.trim_start();

        let in_continuation = nesting_depth > 0 || pending_backslash;
        let should_trace = !in_continuation
            && !trimmed.is_empty()
            && !trimmed.starts_with('#')
            && !trimmed.starts_with('@')
            && !trimmed.starts_with('"')
            && !trimmed.starts_with('\'')
            && !trimmed.starts_with(']')
            && !trimmed.starts_with(')')
            && !trimmed.starts_with('}')
            && !trimmed.starts_with(',')
            && !trimmed.starts_with('+')
            && !trimmed.starts_with('-')
            && !trimmed.starts_with('*')
            && !trimmed.starts_with('/')
            && !trimmed.starts_with('%')
            && !trimmed.starts_with("and ")
            && !trimmed.starts_with("or ")
            && !trimmed.starts_with("if ")
            && !trimmed.starts_with("for ")
            && !trimmed.starts_with("elif ")
            && !trimmed.starts_with("else:")
            && !trimmed.starts_with("except")
            && !trimmed.starts_with("finally:")
            && !trimmed.ends_with(':');

        if should_trace {
            let indent_len = line.len() - trimmed.len();
            let indent = &line[..indent_len];
            out.push(format!("{indent}__bioscript_trace__({line_no})"));
        }
        out.push(line.to_owned());

        pending_backslash = ends_with_unescaped_backslash(line);
        nesting_depth = update_nesting_depth(nesting_depth, line);
    }
    if code.ends_with('\n') {
        out.join("\n") + "\n"
    } else {
        out.join("\n")
    }
}

pub(crate) fn ends_with_unescaped_backslash(line: &str) -> bool {
    let trimmed = line.trim_end();
    if !trimmed.ends_with('\\') {
        return false;
    }

    let slash_count = trimmed.chars().rev().take_while(|ch| *ch == '\\').count();
    slash_count % 2 == 1
}

pub(crate) fn update_nesting_depth(mut depth: usize, line: &str) -> usize {
    let mut chars = line.chars().peekable();
    let mut in_single = false;
    let mut in_double = false;

    while let Some(ch) = chars.next() {
        if in_single {
            if ch == '\\' {
                chars.next();
            } else if ch == '\'' {
                in_single = false;
            }
            continue;
        }

        if in_double {
            if ch == '\\' {
                chars.next();
            } else if ch == '"' {
                in_double = false;
            }
            continue;
        }

        match ch {
            '#' => break,
            '\'' => in_single = true,
            '"' => in_double = true,
            '(' | '[' | '{' => depth += 1,
            ')' | ']' | '}' => depth = depth.saturating_sub(1),
            _ => {}
        }
    }

    depth
}
