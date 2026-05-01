use std::collections::{BTreeSet, HashMap};

use bioscript_core::RuntimeError;

use super::normalize_genotype;

mod scan;

pub(crate) use scan::scan_delimited_variants;

const COMMENT_PREFIXES: [&str; 2] = ["#", "//"];
const RSID_ALIASES: &[&str] = &["rsid", "name", "snp", "marker", "id", "snpid"];
const CHROM_ALIASES: &[&str] = &["chromosome", "chr", "chrom"];
const POSITION_ALIASES: &[&str] = &[
    "position",
    "pos",
    "coordinate",
    "basepairposition",
    "basepair",
];
pub(crate) const GENOTYPE_ALIASES: &[&str] = &[
    "genotype",
    "gt",
    "result",
    "results",
    "result1",
    "call",
    "calls",
    "yourcode",
    "code",
    "genotypevalue",
    "variation",
];
const ALLELE1_ALIASES: &[&str] = &["allele1", "allelea", "allele_a", "allele1top"];
const ALLELE2_ALIASES: &[&str] = &["allele2", "alleleb", "allele_b", "allele2top"];

#[derive(Debug, Clone)]
pub(crate) struct ParsedDelimitedRow {
    pub(crate) rsid: Option<String>,
    pub(crate) chrom: Option<String>,
    pub(crate) position: Option<i64>,
    pub(crate) genotype: String,
}

#[derive(Debug, Clone, Copy)]
pub(crate) enum Delimiter {
    Tab,
    Comma,
    Space,
}

pub(crate) fn detect_delimiter(lines: &[String]) -> Delimiter {
    for line in lines {
        let trimmed = line.trim();
        if trimmed.is_empty()
            || COMMENT_PREFIXES
                .iter()
                .any(|prefix| trimmed.starts_with(prefix))
        {
            continue;
        }
        if line.contains('\t') {
            return Delimiter::Tab;
        }
        if line.contains(',') {
            return Delimiter::Comma;
        }
        if trimmed.split_whitespace().count() > 1 {
            return Delimiter::Space;
        }
    }
    Delimiter::Tab
}

#[allow(dead_code)]
pub(crate) struct RowParser {
    delimiter: Delimiter,
    header: Option<Vec<String>>,
    comment_header: Option<Vec<String>>,
    alias_map: HashMap<&'static str, BTreeSet<&'static str>>,
}

#[allow(dead_code)]
impl RowParser {
    pub(crate) fn new(delimiter: Delimiter) -> Self {
        let mut alias_map = HashMap::new();
        alias_map.insert("rsid", RSID_ALIASES.iter().copied().collect());
        alias_map.insert("chromosome", CHROM_ALIASES.iter().copied().collect());
        alias_map.insert("position", POSITION_ALIASES.iter().copied().collect());
        alias_map.insert("genotype", GENOTYPE_ALIASES.iter().copied().collect());
        alias_map.insert("allele1", ALLELE1_ALIASES.iter().copied().collect());
        alias_map.insert("allele2", ALLELE2_ALIASES.iter().copied().collect());
        Self {
            delimiter,
            header: None,
            comment_header: None,
            alias_map,
        }
    }

    pub(crate) fn consume_line(
        &mut self,
        line: &str,
    ) -> Result<Option<(String, String)>, RuntimeError> {
        Ok(self
            .consume_record(line)?
            .and_then(|row| row.rsid.map(|rsid| (rsid, row.genotype))))
    }

    pub(crate) fn consume_record(
        &mut self,
        line: &str,
    ) -> Result<Option<ParsedDelimitedRow>, RuntimeError> {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            return Ok(None);
        }

        let trimmed = strip_bom(trimmed);
        if let Some(prefix) = COMMENT_PREFIXES
            .iter()
            .find(|prefix| trimmed.starts_with(**prefix))
        {
            let candidate = trimmed.trim_start_matches(prefix).trim();
            if !candidate.is_empty() {
                let fields = self.parse_fields(candidate);
                if self.looks_like_header(&fields) {
                    self.comment_header = Some(fields);
                }
            }
            return Ok(None);
        }

        let fields = self.parse_fields(strip_bom(line));
        if fields.is_empty() {
            return Ok(None);
        }

        if self.header.is_none() {
            if self.looks_like_header(&fields) {
                self.header = Some(fields);
                return Ok(None);
            }
            if let Some(header) = self.comment_header.take() {
                self.header = Some(header);
            } else {
                self.header = Some(self.default_header(fields.len()));
            }
        }

        let header = self.header.as_ref().expect("header initialized");
        let mut row_map = HashMap::new();
        for (idx, value) in fields.into_iter().enumerate() {
            if idx >= header.len() {
                continue;
            }
            row_map.insert(normalize_name(&header[idx]), strip_inline_comment(&value));
        }

        let rsid = self
            .lookup(&row_map, "rsid")
            .filter(|value| !value.is_empty());
        let chrom = self
            .lookup(&row_map, "chromosome")
            .filter(|value| !value.is_empty());
        let position = self
            .lookup(&row_map, "position")
            .and_then(|value| value.parse::<i64>().ok());
        if rsid.is_none() && (chrom.is_none() || position.is_none()) {
            return Ok(None);
        }

        let genotype = if let Some(gt) = self.lookup(&row_map, "genotype") {
            gt
        } else {
            let allele1 = self.lookup(&row_map, "allele1").unwrap_or_default();
            let allele2 = self.lookup(&row_map, "allele2").unwrap_or_default();
            format!("{allele1}{allele2}")
        };

        Ok(Some(ParsedDelimitedRow {
            rsid,
            chrom,
            position,
            genotype: normalize_genotype(&genotype),
        }))
    }

    fn parse_fields(&self, line: &str) -> Vec<String> {
        parse_owned_fields(line, self.delimiter)
    }

    fn looks_like_header(&self, fields: &[String]) -> bool {
        fields.first().is_some_and(|first| {
            self.alias_map
                .get("rsid")
                .is_some_and(|aliases| aliases.contains(normalize_name(first).as_str()))
        })
    }

    fn lookup(&self, row_map: &HashMap<String, String>, key: &str) -> Option<String> {
        let aliases = self.alias_map.get(key)?;
        for alias in aliases {
            let key = normalize_name(alias);
            if let Some(value) = row_map.get(&key)
                && !value.is_empty()
            {
                return Some(value.clone());
            }
        }
        None
    }

    pub(crate) fn default_header(&self, field_count: usize) -> Vec<String> {
        let base = ["rsid", "chromosome", "position", "genotype"];
        if field_count <= base.len() {
            base[..field_count]
                .iter()
                .map(|s| (*s).to_owned())
                .collect()
        } else {
            let mut header: Vec<String> = base.iter().map(|s| (*s).to_owned()).collect();
            for idx in 0..(field_count - header.len()) {
                header.push(format!("extra_{idx}"));
            }
            header
        }
    }
}

pub(crate) fn strip_bom(value: &str) -> &str {
    value.strip_prefix('\u{feff}').unwrap_or(value)
}

pub(crate) fn normalize_name(name: &str) -> String {
    name.trim()
        .to_ascii_lowercase()
        .chars()
        .filter(|ch| !matches!(ch, ' ' | '_' | '-'))
        .collect()
}

pub(crate) fn strip_inline_comment(value: &str) -> String {
    for marker in ["#", "//"] {
        if let Some(idx) = value.find(marker) {
            return value[..idx].trim().to_owned();
        }
    }
    value.trim().to_owned()
}

pub(crate) fn split_csv_line(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut current = String::new();
    let mut in_quotes = false;
    let chars = line.chars().peekable();

    for ch in chars {
        match ch {
            '"' => in_quotes = !in_quotes,
            ',' if !in_quotes => {
                fields.push(current.trim().to_owned());
                current.clear();
            }
            _ => current.push(ch),
        }
    }
    fields.push(current.trim().to_owned());
    fields
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct DelimitedColumnIndexes {
    pub(crate) rsid: Option<usize>,
    pub(crate) chrom: Option<usize>,
    pub(crate) position: Option<usize>,
    pub(crate) genotype: Option<usize>,
    pub(crate) allele1: Option<usize>,
    pub(crate) allele2: Option<usize>,
}

pub(crate) fn parse_streaming_row(
    line: &str,
    delimiter: Delimiter,
    column_indexes: &mut Option<DelimitedColumnIndexes>,
    comment_header: &mut Option<Vec<String>>,
) -> Result<Option<ParsedDelimitedRow>, RuntimeError> {
    let trimmed = line.trim();
    if trimmed.is_empty() {
        return Ok(None);
    }

    let trimmed = strip_bom(trimmed);
    if let Some(prefix) = COMMENT_PREFIXES
        .iter()
        .find(|prefix| trimmed.starts_with(**prefix))
    {
        let candidate = trimmed.trim_start_matches(prefix).trim();
        if !candidate.is_empty() {
            let fields = parse_owned_fields(candidate, delimiter);
            if looks_like_header_fields(&fields) {
                *comment_header = Some(fields);
            }
        }
        return Ok(None);
    }

    let fields = parse_owned_fields(strip_bom(line), delimiter);
    if fields.is_empty() {
        return Ok(None);
    }

    if column_indexes.is_none() {
        if looks_like_header_fields(&fields) {
            *column_indexes = Some(build_column_indexes(&fields));
            return Ok(None);
        }
        if let Some(header) = comment_header.take() {
            *column_indexes = Some(build_column_indexes(&header));
        } else {
            *column_indexes = Some(default_column_indexes(fields.len()));
        }
    }

    let indexes = column_indexes.expect("streaming column indexes initialized");
    let rsid = indexes
        .rsid
        .and_then(|idx| fields.get(idx))
        .map(|value| strip_inline_comment(value).trim().to_owned())
        .filter(|value| !value.is_empty());
    let chrom = indexes
        .chrom
        .and_then(|idx| fields.get(idx))
        .map(|value| strip_inline_comment(value).trim().to_owned())
        .filter(|value| !value.is_empty());
    let position = indexes
        .position
        .and_then(|idx| fields.get(idx))
        .and_then(|value| strip_inline_comment(value).trim().parse::<i64>().ok());
    if rsid.is_none() && (chrom.is_none() || position.is_none()) {
        return Ok(None);
    }

    let genotype = if let Some(idx) = indexes.genotype {
        fields
            .get(idx)
            .map(|value| strip_inline_comment(value))
            .unwrap_or_default()
            .clone()
    } else {
        let allele1 = indexes
            .allele1
            .and_then(|idx| fields.get(idx))
            .map(|value| strip_inline_comment(value))
            .unwrap_or_default();
        let allele2 = indexes
            .allele2
            .and_then(|idx| fields.get(idx))
            .map(|value| strip_inline_comment(value))
            .unwrap_or_default();
        format!("{allele1}{allele2}")
    };

    Ok(Some(ParsedDelimitedRow {
        rsid,
        chrom,
        position,
        genotype: normalize_genotype(&genotype),
    }))
}

fn parse_owned_fields(line: &str, delimiter: Delimiter) -> Vec<String> {
    match delimiter {
        Delimiter::Tab => line
            .split('\t')
            .map(|field| field.trim().to_owned())
            .collect(),
        Delimiter::Space => line.split_whitespace().map(str::to_owned).collect(),
        Delimiter::Comma => split_csv_line(line),
    }
}

pub(crate) fn looks_like_header_fields(fields: &[String]) -> bool {
    fields
        .first()
        .is_some_and(|first| RSID_ALIASES.contains(&normalize_name(first).as_str()))
}

pub(crate) fn build_column_indexes(header: &[String]) -> DelimitedColumnIndexes {
    DelimitedColumnIndexes {
        rsid: find_header_index(header, RSID_ALIASES),
        chrom: find_header_index(header, CHROM_ALIASES),
        position: find_header_index(header, POSITION_ALIASES),
        genotype: find_header_index(header, GENOTYPE_ALIASES),
        allele1: find_header_index(header, ALLELE1_ALIASES),
        allele2: find_header_index(header, ALLELE2_ALIASES),
    }
}

pub(crate) fn default_column_indexes(field_count: usize) -> DelimitedColumnIndexes {
    DelimitedColumnIndexes {
        rsid: (field_count > 0).then_some(0),
        chrom: (field_count > 1).then_some(1),
        position: (field_count > 2).then_some(2),
        genotype: (field_count > 3).then_some(3),
        allele1: None,
        allele2: None,
    }
}

pub(crate) fn find_header_index(header: &[String], aliases: &[&str]) -> Option<usize> {
    header.iter().position(|field| {
        aliases
            .iter()
            .any(|alias| normalize_name(field) == normalize_name(alias))
    })
}
