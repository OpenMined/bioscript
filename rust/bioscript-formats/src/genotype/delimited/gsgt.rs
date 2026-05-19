//! Illumina `GenomeStudio` GSGT Final Report ("Carigenetics") support.
//!
//! Maps GSGT rows into the existing `ParsedDelimitedRow` so the standard
//! rsid-first / `(chrom,pos)` locus matching engine handles the rest. See
//! `illumina.md` for the format spec and the quirks this module encodes.

use bioscript_core::RuntimeError;

use super::{ParsedDelimitedRow, normalize_name, sanitize_evidence_line};
use crate::genotype::normalize_genotype;

/// True when the sampled lines are a GSGT Final Report: the first non-empty
/// line is `[Header]` (case-insensitive).
pub(crate) fn lines_look_like_gsgt(lines: &[String]) -> bool {
    for line in lines {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        return trimmed.eq_ignore_ascii_case("[header]");
    }
    false
}

/// Extract the first `rs\d+` substring from a GSGT `SNP Name` probe id.
///
/// `BOT-rs1135675` -> `rs1135675`, `rs111647200_ilmndup1` -> `rs111647200`,
/// `1:103380393` -> `None` (the embedded coordinate is *not* an rsid).
pub(crate) fn extract_rsid(snp_name: &str) -> Option<String> {
    let bytes = snp_name.as_bytes();
    let mut i = 0;
    while i + 2 < bytes.len() {
        let is_rs = (bytes[i] == b'r' || bytes[i] == b'R')
            && (bytes[i + 1] == b's' || bytes[i + 1] == b'S');
        if is_rs && bytes[i + 2].is_ascii_digit() {
            let mut j = i + 2;
            while j < bytes.len() && bytes[j].is_ascii_digit() {
                j += 1;
            }
            return Some(format!("rs{}", &snp_name[i + 2..j]));
        }
        i += 1;
    }
    None
}

#[derive(Debug, Clone, Copy)]
struct GsgtColumns {
    snp_name: usize,
    chr: usize,
    position: usize,
    allele1_plus: usize,
    allele2_plus: usize,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Phase {
    /// Inside the `[Header]` metadata block, before `[Data]`.
    Meta,
    /// `[Data]` seen; the next non-empty line is the column header.
    ExpectHeader,
    /// Column header parsed; remaining lines are data rows.
    Body,
}

/// Streaming GSGT row parser. Feed it every line of the file in order; it
/// skips the metadata block, learns the columns, then yields one
/// `ParsedDelimitedRow` per usable data row.
#[derive(Debug)]
pub(crate) struct GsgtParser {
    phase: Phase,
    cols: Option<GsgtColumns>,
}

impl GsgtParser {
    pub(crate) fn new() -> Self {
        Self {
            phase: Phase::Meta,
            cols: None,
        }
    }

    pub(crate) fn consume(
        &mut self,
        line: &str,
    ) -> Result<Option<ParsedDelimitedRow>, RuntimeError> {
        let trimmed = line.trim();
        match self.phase {
            Phase::Meta => {
                if trimmed.eq_ignore_ascii_case("[data]") {
                    self.phase = Phase::ExpectHeader;
                }
                Ok(None)
            }
            Phase::ExpectHeader => {
                if trimmed.is_empty() {
                    return Ok(None);
                }
                self.cols = Some(resolve_columns(&split_tab(line))?);
                self.phase = Phase::Body;
                Ok(None)
            }
            Phase::Body => {
                if trimmed.is_empty() {
                    return Ok(None);
                }
                let cols = self.cols.expect("columns resolved before body");
                let fields = split_tab(line);
                Ok(extract_row(line, &fields, &cols))
            }
        }
    }
}

fn split_tab(line: &str) -> Vec<&str> {
    line.trim_end_matches(['\n', '\r'])
        .split('\t')
        .map(str::trim)
        .collect()
}

fn resolve_columns(header: &[&str]) -> Result<GsgtColumns, RuntimeError> {
    let find = |want: &str| header.iter().position(|h| normalize_name(h) == want);
    let snp_name = find("snpname");
    let chr = find("chr").or_else(|| find("chromosome")).or_else(|| find("chrom"));
    let position = find("position").or_else(|| find("pos"));
    let allele1_plus = find("allele1plus");
    let allele2_plus = find("allele2plus");
    match (snp_name, chr, position, allele1_plus, allele2_plus) {
        (Some(snp_name), Some(chr), Some(position), Some(allele1_plus), Some(allele2_plus)) => {
            Ok(GsgtColumns {
                snp_name,
                chr,
                position,
                allele1_plus,
                allele2_plus,
            })
        }
        _ => Err(RuntimeError::Unsupported(
            "GSGT Final Report missing required columns (SNP Name, Chr, Position, Allele1 - Plus, Allele2 - Plus)"
                .to_owned(),
        )),
    }
}

fn extract_row(
    line: &str,
    fields: &[&str],
    cols: &GsgtColumns,
) -> Option<ParsedDelimitedRow> {
    let snp_name = fields.get(cols.snp_name).copied().unwrap_or("");
    let rsid = extract_rsid(snp_name);

    let chrom_raw = fields.get(cols.chr).copied().unwrap_or("").trim();
    // Skip unplaced markers (Chr/Position == 0).
    if chrom_raw.is_empty() || chrom_raw == "0" {
        return None;
    }
    let position = fields
        .get(cols.position)
        .and_then(|value| value.trim().parse::<i64>().ok())
        .filter(|pos| *pos != 0)?;

    if rsid.is_none() {
        // No rsid and we still have a locus; keep it for the locus fallback.
        // (chrom/position validated above.)
    }

    let a1 = fields.get(cols.allele1_plus).copied().unwrap_or("").trim();
    let a2 = fields.get(cols.allele2_plus).copied().unwrap_or("").trim();
    // GSGT no-call is a single `-` per Plus allele column.
    let genotype = if (a1 == "-" || a1.is_empty()) && (a2 == "-" || a2.is_empty()) {
        "--".to_owned()
    } else {
        normalize_genotype(&format!("{a1}{a2}"))
    };

    Some(ParsedDelimitedRow {
        rsid,
        chrom: Some(chrom_raw.to_owned()),
        position: Some(position),
        genotype,
        raw_line: sanitize_evidence_line(line),
    })
}

/// Whether a normalized genotype string represents a no-call.
pub(crate) fn is_no_call(genotype: &str) -> bool {
    genotype == "--" || genotype.is_empty()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lines_look_like_gsgt_detects_header_marker() {
        assert!(lines_look_like_gsgt(&[
            String::new(),
            "[Header]".to_owned()
        ]));
        assert!(lines_look_like_gsgt(&["[HEADER]".to_owned()]));
        assert!(!lines_look_like_gsgt(&[
            "# Dynamic DNA".to_owned(),
            "[Header]".to_owned()
        ]));
        assert!(!lines_look_like_gsgt(&[
            "rsid\tchromosome\tposition\tgenotype".to_owned()
        ]));
    }

    #[test]
    fn extract_rsid_recovers_rs_from_every_probe_form() {
        assert_eq!(extract_rsid("rs11466023").as_deref(), Some("rs11466023"));
        assert_eq!(extract_rsid("BOT-rs1135675").as_deref(), Some("rs1135675"));
        assert_eq!(
            extract_rsid("rs111647200_ilmndup1").as_deref(),
            Some("rs111647200")
        );
        assert_eq!(extract_rsid("GSA-rs61660502").as_deref(), Some("rs61660502"));
        assert_eq!(extract_rsid("seq-rs786202193").as_deref(), Some("rs786202193"));
        // No rs id: chr:pos, CNV, MNV, vendor/clinical.
        assert_eq!(extract_rsid("1:103380393"), None);
        assert_eq!(extract_rsid("1:159174749-C-T"), None);
        assert_eq!(extract_rsid("1:110228436_CNV_GSTM1"), None);
        assert_eq!(extract_rsid("1:45332290_MNV"), None);
        assert_eq!(extract_rsid("DICER1-chr14-95596479"), None);
        assert_eq!(extract_rsid("GALC:c.2002A>C"), None);
    }

    fn body_parser() -> GsgtParser {
        let mut p = GsgtParser::new();
        for line in [
            "[Header]",
            "GSGT Version\t2.0.5",
            "[Data]",
            "Sample ID\tSample Name\tSNP Name\tSNP\tChr\tPosition\tAllele1 - Top\tAllele2 - Top\tAllele1 - Plus\tAllele2 - Plus\tPlus/Minus Strand",
        ] {
            assert!(p.consume(line).unwrap().is_none());
        }
        p
    }

    fn row(p: &mut GsgtParser, snp_name: &str, snp: &str, chr: &str, pos: &str, a1: &str, a2: &str) -> Option<ParsedDelimitedRow> {
        let line = format!("S1\t\t{snp_name}\t{snp}\t{chr}\t{pos}\tX\tX\t{a1}\t{a2}\t+");
        p.consume(&line).unwrap()
    }

    #[test]
    fn parses_plus_columns_and_snp_name_not_design_column() {
        let mut p = body_parser();
        let r = row(&mut p, "rs9000001", "[A/G]", "1", "1000", "G", "G").unwrap();
        assert_eq!(r.rsid.as_deref(), Some("rs9000001"));
        assert_eq!(r.chrom.as_deref(), Some("1"));
        assert_eq!(r.position, Some(1000));
        assert_eq!(r.genotype, "GG");
    }

    #[test]
    fn design_snp_column_is_never_used_as_rsid() {
        let mut p = body_parser();
        // SNP Name has no rs id; SNP design column is [A/G]. rsid must be None.
        let r = row(&mut p, "1:7000_CNV_GENEX", "[A/G]", "1", "7000", "A", "A").unwrap();
        assert_eq!(r.rsid, None);
        assert_eq!(r.position, Some(7000));
        assert_eq!(r.genotype, "AA");
    }

    #[test]
    fn no_call_dash_maps_to_double_dash_and_indels_pass_through() {
        let mut p = body_parser();
        let nc = row(&mut p, "rs9000011", "[A/G]", "4", "9000", "-", "-").unwrap();
        assert_eq!(nc.genotype, "--");
        let indel = row(&mut p, "rs9000012", "[I/D]", "4", "9100", "I", "D").unwrap();
        assert!(indel.genotype != "--");
    }

    #[test]
    fn skips_unplaced_chr_or_position_zero() {
        let mut p = body_parser();
        assert!(row(&mut p, "rs9000009", "[A/G]", "0", "9000", "A", "G").is_none());
        assert!(row(&mut p, "rs9000010", "[A/G]", "4", "0", "A", "G").is_none());
    }

    #[test]
    fn empty_rsid_kept_when_locus_present() {
        let mut p = body_parser();
        let r = row(&mut p, "2:6000", "[C/T]", "2", "6000", "C", "T").unwrap();
        assert_eq!(r.rsid, None);
        assert_eq!(r.chrom.as_deref(), Some("2"));
        assert_eq!(r.position, Some(6000));
        assert_eq!(r.genotype, "CT");
    }
}
