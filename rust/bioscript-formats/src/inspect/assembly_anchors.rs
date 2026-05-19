//! rsID / locus anchor genome-build detector.
//!
//! Used as the fallback when a genotype/VCF file declares no build in its
//! metadata (header build line, VCF `##reference`/contig lengths). Compares
//! `Chr`+`Position` of well-known SNPs against a table of their GRCh36 /
//! GRCh37 / GRCh38 coordinates; the build with the dominant exact-position
//! vote wins. Positions discriminate the build; alleles only corroborate a
//! locus hit (strand-agnostic) so a file with no usable rsIDs (e.g. a VCF
//! with `ID=.`, or GSGT non-rs probes) can still be placed.
//!
//! The anchor table was self-sourced from declared-build ground-truth files
//! (GRCh36 = 23andMe v2/v3, GRCh37 = AncestryDNA + 23andMe v4/v5,
//! GRCh38 = Dynamic DNA / GSGT) and validated to classify every cached
//! vendor file correctly; see `illumina.md`. GRCh36 is present only so a
//! build-36 file does not get misread as 37/38 — `decide` only ever returns
//! GRCh37, GRCh38, or `None` (honest unknown).

#![allow(clippy::doc_markdown)]

use bioscript_core::Assembly;

/// `(rsid, chrom, grch36_pos, grch37_pos, grch38_pos, plus_strand_alleles)`.
type Anchor = (&'static str, &'static str, i64, i64, i64, &'static str);

const ANCHORS: &[Anchor] = &[
    ("rs1801133", "1", 11_778_965, 11_856_378, 11_796_321, "AG"),
    ("rs6025", "1", 167_785_673, 169_519_049, 169_549_811, "C"),
    ("rs3131972", "1", 742_584, 752_721, 817_341, "AG"),
    ("rs4988235", "2", 136_325_116, 136_608_646, 135_851_076, "G"),
    ("rs1801282", "3", 12_368_125, 12_393_125, 12_351_626, "C"),
    ("rs53576", "3", 8_779_371, 8_804_371, 8_762_685, "AG"),
    ("rs1229984", "4", 100_458_342, 100_239_319, 99_318_162, "C"),
    ("rs1801260", "4", 55_996_126, 56_301_369, 55_435_202, "AG"),
    ("rs16891982", "5", 33_987_450, 33_951_693, 33_951_588, "CG"),
    ("rs1042713", "5", 148_186_633, 148_206_440, 148_826_877, "G"),
    ("rs1799971", "6", 154_402_490, 154_360_797, 154_039_662, "A"),
    ("rs662", "7", 94_775_382, 94_937_446, 95_308_134, "CT"),
    ("rs505922", "9", 135_139_050, 136_149_229, 133_273_813, "CT"),
    (
        "rs7903146",
        "10",
        114_748_339,
        114_758_349,
        112_998_590,
        "C",
    ),
    (
        "rs1800497",
        "11",
        112_776_038,
        113_270_828,
        113_400_106,
        "G",
    ),
    ("rs671", "12", 110_726_149, 112_241_766, 111_803_962, "G"),
    ("rs12913832", "15", 26_039_213, 28_365_618, 28_120_472, "AG"),
    ("rs1051730", "15", 76_681_394, 78_894_339, 78_601_997, "G"),
    ("rs16969968", "15", 76_669_980, 78_882_925, 78_590_583, "G"),
    ("rs762551", "15", 72_828_970, 75_041_917, 74_749_576, "AC"),
    ("rs1042522", "17", 7_520_197, 7_579_472, 7_676_154, "CG"),
    ("rs429358", "19", 50_103_781, 45_411_941, 44_908_684, "CT"),
    ("rs7412", "19", 50_103_919, 45_412_079, 44_908_822, "C"),
    ("rs4680", "22", 18_331_271, 19_951_271, 19_963_748, "AG"),
    ("rs1800896", "1", 205_013_520, 206_946_897, 206_773_552, "T"),
    ("rs1800795", "7", 22_733_170, 22_766_645, 22_727_026, "G"),
    ("rs5186", "3", 149_942_678, 148_459_988, 148_742_201, "A"),
];

const RSID_WEIGHT: u32 = 3;
const LOCUS_WEIGHT: u32 = 1;
const MIN_ANCHORS: u32 = 4;
const MIN_DOMINANCE: f64 = 0.90;

fn complement(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(b'T'),
        b'T' => Some(b'A'),
        b'C' => Some(b'G'),
        b'G' => Some(b'C'),
        _ => None,
    }
}

/// Normalize a chromosome label to bare form (`chr7` -> `7`, `23` -> `X`).
fn norm_chrom(raw: &str) -> String {
    let c = raw.trim();
    let c = c
        .strip_prefix("chr")
        .or_else(|| c.strip_prefix("CHR"))
        .unwrap_or(c);
    match c.to_ascii_uppercase().as_str() {
        "23" => "X".to_owned(),
        "24" => "Y".to_owned(),
        "25" => "XY".to_owned(),
        "26" | "M" => "MT".to_owned(),
        _ => c.to_ascii_uppercase(),
    }
}

/// Observed ACGT letters in a genotype must be a subset of the anchor's
/// allele set or its complement (strand-agnostic). Empty/no-call passes.
fn allele_ok(genotype: &str, anchor_alleles: &str) -> bool {
    let observed: Vec<u8> = genotype
        .bytes()
        .map(|b| b.to_ascii_uppercase())
        .filter(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
        .collect();
    if observed.is_empty() {
        return true;
    }
    let fwd: Vec<u8> = anchor_alleles
        .bytes()
        .map(|b| b.to_ascii_uppercase())
        .collect();
    let rev: Vec<u8> = fwd.iter().filter_map(|b| complement(*b)).collect();
    let subset = |set: &[u8]| observed.iter().all(|o| set.contains(o));
    subset(&fwd) || subset(&rev)
}

/// Accumulates anchor votes from a streamed genotype/VCF file.
#[derive(Debug, Default)]
pub(crate) struct AssemblyAnchorScorer {
    rs: [u32; 3],
    loc: [u32; 3],
}

impl AssemblyAnchorScorer {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// Observe one parsed row. `rsid` may be empty (locus-only vote).
    pub(crate) fn observe(&mut self, rsid: &str, chrom: &str, pos: i64, genotype: &str) {
        let chrom = norm_chrom(chrom);
        for (a_rsid, a_chrom, p36, p37, p38, alleles) in ANCHORS {
            if chrom != *a_chrom {
                continue;
            }
            let build = if pos == *p36 {
                0
            } else if pos == *p37 {
                1
            } else if pos == *p38 {
                2
            } else {
                continue;
            };
            if rsid == *a_rsid {
                self.rs[build] += 1;
            } else if allele_ok(genotype, alleles) {
                self.loc[build] += 1;
            }
            return;
        }
    }

    fn scores(&self) -> [u32; 3] {
        [
            self.rs[0] * RSID_WEIGHT + self.loc[0] * LOCUS_WEIGHT,
            self.rs[1] * RSID_WEIGHT + self.loc[1] * LOCUS_WEIGHT,
            self.rs[2] * RSID_WEIGHT + self.loc[2] * LOCUS_WEIGHT,
        ]
    }

    /// Resolve to a supported assembly, or `None` when there is too little
    /// signal, the vote is split, or a build-36 file dominates (unsupported).
    pub(crate) fn decide(&self) -> Option<Assembly> {
        let anchors: u32 = self.rs.iter().chain(self.loc.iter()).sum();
        if anchors < MIN_ANCHORS {
            return None;
        }
        let s = self.scores();
        let total: u32 = s.iter().sum();
        if total == 0 {
            return None;
        }
        let (best_idx, &best) = s
            .iter()
            .enumerate()
            .max_by_key(|(_, v)| **v)
            .expect("three scores");
        if f64::from(best) / f64::from(total) < MIN_DOMINANCE {
            return None;
        }
        match best_idx {
            1 => Some(Assembly::Grch37),
            2 => Some(Assembly::Grch38),
            // best_idx == 0 -> GRCh36, which we do not model: honest unknown.
            _ => None,
        }
    }
}

/// Vote an assembly from an in-memory text buffer (genotype text or VCF),
/// vendor-agnostic. Handles GSGT `[Header]/[Data]`, VCF, and flat
/// rsid/chrom/pos/genotype delimited rows. Returns `None` unless the anchor
/// vote is confident — used only as the metadata-absent fallback.
#[allow(clippy::items_after_statements, clippy::many_single_char_names)]
pub(crate) fn assembly_from_text_bytes(bytes: &[u8]) -> Option<Assembly> {
    enum Mode {
        Auto,
        GsgtMeta,
        GsgtHdr,
        GsgtBody(usize, usize, usize, usize, usize),
        Vcf,
    }
    let text = String::from_utf8_lossy(bytes);
    let mut scorer = AssemblyAnchorScorer::new();
    let mut mode = Mode::Auto;
    for raw in text.lines() {
        let line = raw.trim_end_matches('\r');
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        match mode {
            Mode::Auto if t.eq_ignore_ascii_case("[header]") => {
                mode = Mode::GsgtMeta;
            }
            Mode::GsgtMeta => {
                if t.eq_ignore_ascii_case("[data]") {
                    mode = Mode::GsgtHdr;
                }
            }
            Mode::GsgtHdr => {
                let norm = |s: &str| s.trim().to_ascii_lowercase().replace([' ', '-', '_'], "");
                let h: Vec<String> = line.split('\t').map(norm).collect();
                let idx = |w: &str| h.iter().position(|x| x == w);
                match (
                    idx("snpname"),
                    idx("chr"),
                    idx("position"),
                    idx("allele1plus"),
                    idx("allele2plus"),
                ) {
                    (Some(s), Some(c), Some(p), Some(a1), Some(a2)) => {
                        mode = Mode::GsgtBody(s, c, p, a1, a2);
                    }
                    _ => return None,
                }
            }
            Mode::GsgtBody(s, c, p, a1, a2) => {
                let f: Vec<&str> = line.split('\t').collect();
                if f.len() <= a2.max(p).max(c).max(s) {
                    continue;
                }
                if let Ok(pos) = f[p].trim().parse::<i64>() {
                    let rsid = extract_rs(f[s]);
                    scorer.observe(
                        rsid.as_deref().unwrap_or(""),
                        f[c].trim(),
                        pos,
                        &format!("{}{}", f[a1].trim(), f[a2].trim()),
                    );
                }
            }
            Mode::Vcf => {
                let f: Vec<&str> = t.split('\t').collect();
                if f.len() >= 5
                    && let Ok(pos) = f[1].trim().parse::<i64>()
                {
                    let rsid = if f[2].starts_with("rs") { f[2] } else { "" };
                    scorer.observe(rsid, f[0], pos, &format!("{}{}", f[3], f[4]));
                }
            }
            Mode::Auto => {
                if t.starts_with("##") {
                    continue;
                }
                if t.starts_with("#CHROM\t") {
                    mode = Mode::Vcf;
                    continue;
                }
                if t.starts_with('#') || t.starts_with("//") {
                    continue;
                }
                let f: Vec<&str> = if line.contains('\t') {
                    line.split('\t').collect()
                } else {
                    line.split(',').collect()
                };
                if f.len() < 3 {
                    continue;
                }
                let c0 = f[0].trim().trim_matches('"').to_ascii_lowercase();
                if matches!(c0.as_str(), "rsid" | "rs id" | "snp" | "name" | "snpname") {
                    continue; // column header row
                }
                if let Ok(pos) = f[2].trim().trim_matches('"').parse::<i64>() {
                    let id = f[0].trim().trim_matches('"');
                    let rsid = if id.starts_with("rs") || id.starts_with("RS") {
                        id
                    } else {
                        ""
                    };
                    let gt: String = f[3..].iter().map(|x| x.trim().trim_matches('"')).collect();
                    scorer.observe(rsid, f[1].trim().trim_matches('"'), pos, &gt);
                }
            }
        }
    }
    scorer.decide()
}

fn extract_rs(name: &str) -> Option<String> {
    let b = name.as_bytes();
    let mut i = 0;
    while i + 2 < b.len() {
        if (b[i] | 32) == b'r' && (b[i + 1] | 32) == b's' && b[i + 2].is_ascii_digit() {
            let mut j = i + 2;
            while j < b.len() && b[j].is_ascii_digit() {
                j += 1;
            }
            return Some(format!("rs{}", &name[i + 2..j]));
        }
        i += 1;
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn feed(scorer: &mut AssemblyAnchorScorer, build: usize, with_rsid: bool) {
        for (rsid, chrom, p36, p37, p38, alleles) in ANCHORS {
            let pos = [*p36, *p37, *p38][build];
            let id = if with_rsid { *rsid } else { "" };
            scorer.observe(id, chrom, pos, alleles);
        }
    }

    #[test]
    fn rsid_votes_resolve_grch37_and_grch38() {
        let mut g37 = AssemblyAnchorScorer::new();
        feed(&mut g37, 1, true);
        assert_eq!(g37.decide(), Some(Assembly::Grch37));

        let mut g38 = AssemblyAnchorScorer::new();
        feed(&mut g38, 2, true);
        assert_eq!(g38.decide(), Some(Assembly::Grch38));
    }

    #[test]
    fn locus_only_no_rsid_still_resolves() {
        // VCF-style: no usable rsIDs, alleles corroborate the position.
        let mut g38 = AssemblyAnchorScorer::new();
        feed(&mut g38, 2, false);
        assert_eq!(g38.decide(), Some(Assembly::Grch38));
    }

    #[test]
    fn build36_file_is_honest_unknown_not_misread() {
        let mut g36 = AssemblyAnchorScorer::new();
        feed(&mut g36, 0, true);
        assert_eq!(g36.decide(), None);
    }

    #[test]
    fn too_few_anchors_is_unknown() {
        let mut s = AssemblyAnchorScorer::new();
        s.observe("rs429358", "19", 45_411_941, "CT");
        s.observe("rs7412", "19", 45_412_079, "C");
        assert_eq!(s.decide(), None);
    }

    #[test]
    fn split_vote_is_unknown_not_a_guess() {
        let mut s = AssemblyAnchorScorer::new();
        for (i, (rsid, chrom, p36, p37, p38, al)) in ANCHORS.iter().enumerate() {
            let pos = if i % 2 == 0 { *p37 } else { *p38 };
            let _ = (p36,);
            s.observe(rsid, chrom, pos, al);
        }
        assert_eq!(s.decide(), None);
    }

    #[test]
    fn chrom_prefix_and_strand_complement_are_handled() {
        let mut s = AssemblyAnchorScorer::new();
        // "chr19" prefix + complement-strand alleles for rs429358 (CT -> GA).
        for (rsid, chrom, _p36, p37, _p38, _al) in ANCHORS {
            s.observe(rsid, &format!("chr{chrom}"), *p37, "GA");
        }
        assert_eq!(s.decide(), Some(Assembly::Grch37));
    }
}
