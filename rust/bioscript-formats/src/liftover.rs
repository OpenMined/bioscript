use std::{
    cmp::Ordering,
    collections::HashMap,
    fs::File,
    io::{self, BufRead, BufReader, Cursor, Write},
    path::Path,
};

use bioscript_core::{Assembly, RuntimeError};
use flate2::read::GzDecoder;

const HG19_TO_HG38_CHAIN_GZ: &[u8] = include_bytes!("../assets/liftover/hg19ToHg38.over.chain.gz");

const PRIMARY_CHROMS: [&str; 25] = [
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "MT",
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LiftOverOptions {
    pub primary_only: bool,
}

impl Default for LiftOverOptions {
    fn default() -> Self {
        Self { primary_only: true }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LiftedLocus {
    pub chrom: String,
    pub pos: i64,
    pub strand: char,
    pub score: i64,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct LiftoverStats {
    pub total_markers: u64,
    pub mapped: u64,
    pub unmapped: u64,
    pub reverse_strand_genotypes: u64,
}

#[derive(Debug, Clone)]
pub struct ChainIndex {
    from: Assembly,
    to: Assembly,
    blocks_by_chrom: HashMap<String, ChromBlocks>,
}

#[derive(Debug, Clone)]
struct ChromBlocks {
    starts: Vec<i64>,
    blocks: Vec<ChainBlock>,
}

#[derive(Debug, Clone)]
struct ChainBlock {
    source_start: i64,
    source_end: i64,
    target_start: i64,
    target_name: String,
    target_size: i64,
    target_strand: char,
    score: i64,
    index: usize,
}

#[derive(Debug, Clone)]
struct ChainHeader {
    score: i64,
    source_name: String,
    source_pos: i64,
    target_name: String,
    target_size: i64,
    target_strand: char,
    target_pos: i64,
    index: usize,
}

impl ChainIndex {
    pub fn from_reader<R: BufRead>(
        from: Assembly,
        to: Assembly,
        reader: R,
        options: LiftOverOptions,
    ) -> Result<Self, RuntimeError> {
        if (from, to) != (Assembly::Grch37, Assembly::Grch38) {
            return Err(RuntimeError::Unsupported(
                "liftover currently supports only GRCh37 to GRCh38".to_owned(),
            ));
        }

        let mut blocks_by_chrom: HashMap<String, Vec<ChainBlock>> = HashMap::new();
        let mut current: Option<ChainHeader> = None;
        let mut chain_index = 0usize;

        for raw in reader.lines() {
            let raw =
                raw.map_err(|err| RuntimeError::Io(format!("failed to read chain: {err}")))?;
            let line = raw.trim();
            if line.is_empty() {
                current = None;
                continue;
            }

            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.first() == Some(&"chain") {
                current = parse_chain_header(&parts, chain_index, options.primary_only)?;
                if current.is_some() {
                    chain_index += 1;
                }
                continue;
            }

            let Some(header) = current.as_mut() else {
                continue;
            };
            let size = parse_i64(parts[0], "chain block size")?;
            let source_start = header.source_pos;
            let target_start = header.target_pos;
            let source_end = source_start + size;

            blocks_by_chrom
                .entry(header.source_name.clone())
                .or_default()
                .push(ChainBlock {
                    source_start,
                    source_end,
                    target_start,
                    target_name: header.target_name.clone(),
                    target_size: header.target_size,
                    target_strand: header.target_strand,
                    score: header.score,
                    index: header.index,
                });

            if parts.len() == 3 {
                header.source_pos = source_end + parse_i64(parts[1], "chain dt")?;
                header.target_pos = target_start + size + parse_i64(parts[2], "chain dq")?;
            } else {
                current = None;
            }
        }

        let blocks_by_chrom = blocks_by_chrom
            .into_iter()
            .map(|(chrom, mut blocks)| {
                blocks.sort_by(|a, b| {
                    a.source_start
                        .cmp(&b.source_start)
                        .then_with(|| b.score.cmp(&a.score))
                        .then_with(|| a.index.cmp(&b.index))
                });
                let starts = blocks.iter().map(|block| block.source_start).collect();
                (chrom, ChromBlocks { starts, blocks })
            })
            .collect();

        Ok(Self {
            from,
            to,
            blocks_by_chrom,
        })
    }

    pub fn lookup(
        &self,
        from: Assembly,
        to: Assembly,
        chrom: &str,
        pos1: i64,
    ) -> Option<LiftedLocus> {
        if from != self.from || to != self.to || pos1 <= 0 {
            return None;
        }

        let source_chrom = chain_chrom(chrom);
        let chrom_blocks = self
            .blocks_by_chrom
            .get(&source_chrom)
            .or_else(|| self.blocks_by_chrom.get(chrom))?;
        let pos0 = pos1 - 1;
        let mut i = chrom_blocks.starts.partition_point(|start| *start <= pos0);
        if i == 0 {
            return None;
        }
        i -= 1;

        let mut best: Option<&ChainBlock> = None;
        loop {
            let block = &chrom_blocks.blocks[i];
            if block.source_start <= pos0 && pos0 < block.source_end {
                if best.is_none_or(|candidate| compare_blocks(block, candidate).is_gt()) {
                    best = Some(block);
                }
            }

            if let Some(candidate) = best
                && chrom_blocks.starts[i] < candidate.source_start - 1_000_000
            {
                break;
            }
            if best.is_none() && pos0 - chrom_blocks.starts[i] > 10_000_000 {
                break;
            }
            if i == 0 {
                break;
            }
            i -= 1;
        }

        let best = best?;
        let offset = pos0 - best.source_start;
        let target_pos0 = if best.target_strand == '+' {
            best.target_start + offset
        } else {
            best.target_size - (best.target_start + offset) - 1
        };

        Some(LiftedLocus {
            chrom: canon_chrom(&best.target_name),
            pos: target_pos0 + 1,
            strand: best.target_strand,
            score: best.score,
        })
    }
}

pub fn grch37_to_grch38_chain() -> Result<ChainIndex, RuntimeError> {
    let reader = BufReader::new(GzDecoder::new(Cursor::new(HG19_TO_HG38_CHAIN_GZ)));
    ChainIndex::from_reader(
        Assembly::Grch37,
        Assembly::Grch38,
        reader,
        LiftOverOptions::default(),
    )
}

pub fn convert_23andme_grch37_to_grch38(
    input_path: &Path,
    output_path: &Path,
    unmapped_path: &Path,
) -> Result<LiftoverStats, RuntimeError> {
    let chain = grch37_to_grch38_chain()?;
    convert_23andme_with_chain(input_path, output_path, unmapped_path, &chain)
}

pub fn convert_23andme_with_chain(
    input_path: &Path,
    output_path: &Path,
    unmapped_path: &Path,
    chain: &ChainIndex,
) -> Result<LiftoverStats, RuntimeError> {
    let input = File::open(input_path).map_err(|err| {
        RuntimeError::Io(format!("failed to open {}: {err}", input_path.display()))
    })?;
    let mut output = File::create(output_path).map_err(|err| {
        RuntimeError::Io(format!("failed to create {}: {err}", output_path.display()))
    })?;
    let mut unmapped = File::create(unmapped_path).map_err(|err| {
        RuntimeError::Io(format!(
            "failed to create unmapped report {}: {err}",
            unmapped_path.display()
        ))
    })?;

    convert_23andme_reader_with_chain(BufReader::new(input), &mut output, &mut unmapped, chain)
}

pub fn convert_23andme_reader_with_chain<R: BufRead, W: Write, U: Write>(
    reader: R,
    output: &mut W,
    unmapped: &mut U,
    chain: &ChainIndex,
) -> Result<LiftoverStats, RuntimeError> {
    let mut stats = LiftoverStats::default();
    let mut added_note = false;
    writeln!(unmapped, "rsid\tchromosome\tposition\tgenotype\treason")
        .map_err(write_error("unmapped report"))?;

    for raw in reader.lines() {
        let raw =
            raw.map_err(|err| RuntimeError::Io(format!("failed to read 23andMe txt: {err}")))?;
        let line = raw.trim_end_matches(['\n', '\r']);
        if line.starts_with('#') {
            if !added_note && line.contains("reference human assembly build 37") {
                writeln!(
                    output,
                    "# Coordinates lifted from reference human assembly build 37 to build 38 / GRCh38."
                )
                .map_err(write_error("lifted 23andMe output"))?;
                writeln!(
                    output,
                    "# Genotype bases were reverse-complemented for markers mapping through reverse-strand chains."
                )
                .map_err(write_error("lifted 23andMe output"))?;
                added_note = true;
            }
            writeln!(
                output,
                "{}",
                line.replace("build 37", "build 38 / GRCh38")
                    .replace("Annotation Release 104", "GRCh38")
            )
            .map_err(write_error("lifted 23andMe output"))?;
            continue;
        }

        if line.is_empty() {
            writeln!(output).map_err(write_error("lifted 23andMe output"))?;
            continue;
        }

        stats.total_markers += 1;
        let mut fields: Vec<String> = line.split('\t').map(ToOwned::to_owned).collect();
        if fields.len() < 4 {
            stats.unmapped += 1;
            writeln!(unmapped, "{line}\tmalformed").map_err(write_error("unmapped report"))?;
            continue;
        }

        let pos = match fields[2].parse::<i64>() {
            Ok(pos) => pos,
            Err(_) => {
                stats.unmapped += 1;
                writeln!(
                    unmapped,
                    "{}\t{}\t{}\t{}\tnon_integer_position",
                    fields[0], fields[1], fields[2], fields[3]
                )
                .map_err(write_error("unmapped report"))?;
                continue;
            }
        };

        let Some(lifted) = chain.lookup(Assembly::Grch37, Assembly::Grch38, &fields[1], pos) else {
            stats.unmapped += 1;
            writeln!(
                unmapped,
                "{}\t{}\t{}\t{}\tno_primary_mapping",
                fields[0], fields[1], fields[2], fields[3]
            )
            .map_err(write_error("unmapped report"))?;
            continue;
        };

        if lifted.strand == '-' {
            fields[3] = revcomp_genotype(&fields[3]);
            stats.reverse_strand_genotypes += 1;
        }
        fields[1] = lifted.chrom;
        fields[2] = lifted.pos.to_string();
        writeln!(output, "{}", fields.join("\t")).map_err(write_error("lifted 23andMe output"))?;
        stats.mapped += 1;
    }

    Ok(stats)
}

fn parse_chain_header(
    parts: &[&str],
    index: usize,
    primary_only: bool,
) -> Result<Option<ChainHeader>, RuntimeError> {
    if parts.len() < 13 {
        return Err(RuntimeError::Io("malformed chain header".to_owned()));
    }
    let source_name = parts[2].to_owned();
    let source_strand = parts[4];
    let source_pos = parse_i64(parts[5], "chain tStart")?;
    let target_name = parts[7].to_owned();
    let target_size = parse_i64(parts[8], "chain qSize")?;
    let target_strand = parts[9].chars().next().unwrap_or('+');
    let target_pos = parse_i64(parts[10], "chain qStart")?;

    if source_strand != "+" {
        return Ok(None);
    }
    if primary_only
        && (!is_primary_chrom(&canon_chrom(&source_name))
            || !is_primary_chrom(&canon_chrom(&target_name)))
    {
        return Ok(None);
    }

    Ok(Some(ChainHeader {
        score: parse_i64(parts[1], "chain score")?,
        source_name,
        source_pos,
        target_name,
        target_size,
        target_strand,
        target_pos,
        index,
    }))
}

fn compare_blocks(a: &ChainBlock, b: &ChainBlock) -> Ordering {
    a.score.cmp(&b.score).then_with(|| b.index.cmp(&a.index))
}

fn parse_i64(value: &str, label: &str) -> Result<i64, RuntimeError> {
    value
        .parse::<i64>()
        .map_err(|err| RuntimeError::Io(format!("failed to parse {label} '{value}': {err}")))
}

fn canon_chrom(name: &str) -> String {
    let stripped = name.strip_prefix("chr").unwrap_or(name);
    if stripped == "M" {
        "MT".to_owned()
    } else {
        stripped.to_owned()
    }
}

fn chain_chrom(name: &str) -> String {
    let canon = canon_chrom(name);
    if canon == "MT" {
        "chrM".to_owned()
    } else {
        format!("chr{canon}")
    }
}

fn is_primary_chrom(chrom: &str) -> bool {
    PRIMARY_CHROMS.contains(&chrom)
}

fn revcomp_genotype(genotype: &str) -> String {
    genotype
        .chars()
        .map(|ch| match ch {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            'a' => 't',
            'c' => 'g',
            'g' => 'c',
            't' => 'a',
            other => other,
        })
        .collect()
}

fn write_error(label: &'static str) -> impl FnOnce(io::Error) -> RuntimeError {
    move |err| RuntimeError::Io(format!("failed to write {label}: {err}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tiny_chain() -> ChainIndex {
        let chain = "\
chain 100 chr1 1000 + 9 30 chr1 1000 + 99 120 1\n\
10 5 5\n\
10\n\
\n\
chain 90 chr2 1000 + 49 60 chr2 1000 - 100 111 2\n\
10\n\
";
        ChainIndex::from_reader(
            Assembly::Grch37,
            Assembly::Grch38,
            BufReader::new(chain.as_bytes()),
            LiftOverOptions::default(),
        )
        .unwrap()
    }

    #[test]
    fn bundled_chain_maps_known_23andme_marker() {
        let chain = grch37_to_grch38_chain().unwrap();
        let lifted = chain
            .lookup(Assembly::Grch37, Assembly::Grch38, "19", 46_181_392)
            .unwrap();
        assert_eq!(lifted.chrom, "19");
        assert_eq!(lifted.pos, 45_678_134);
        assert_eq!(lifted.strand, '+');
    }

    #[test]
    fn chain_lookup_handles_forward_and_reverse_blocks() {
        let chain = tiny_chain();
        let forward = chain
            .lookup(Assembly::Grch37, Assembly::Grch38, "1", 11)
            .unwrap();
        assert_eq!(forward.chrom, "1");
        assert_eq!(forward.pos, 101);
        assert_eq!(forward.strand, '+');

        let reverse = chain
            .lookup(Assembly::Grch37, Assembly::Grch38, "2", 50)
            .unwrap();
        assert_eq!(reverse.chrom, "2");
        assert_eq!(reverse.pos, 900);
        assert_eq!(reverse.strand, '-');
    }

    #[test]
    fn convert_23andme_lifts_positions_and_reverse_complements() {
        let chain = tiny_chain();
        let input = b"# We are using reference human assembly build 37\n\
rsForward\t1\t11\tAG\n\
rsReverse\t2\t50\tAC\n\
rsMissing\t3\t10\tTT\n";
        let mut out = Vec::new();
        let mut unmapped = Vec::new();
        let stats = convert_23andme_reader_with_chain(
            BufReader::new(&input[..]),
            &mut out,
            &mut unmapped,
            &chain,
        )
        .unwrap();

        assert_eq!(stats.total_markers, 3);
        assert_eq!(stats.mapped, 2);
        assert_eq!(stats.unmapped, 1);
        assert_eq!(stats.reverse_strand_genotypes, 1);

        let out = String::from_utf8(out).unwrap();
        assert!(out.contains("# Coordinates lifted"));
        assert!(out.contains("rsForward\t1\t101\tAG"));
        assert!(out.contains("rsReverse\t2\t900\tTG"));

        let unmapped = String::from_utf8(unmapped).unwrap();
        assert!(unmapped.contains("rsMissing\t3\t10\tTT\tno_primary_mapping"));
    }
}
