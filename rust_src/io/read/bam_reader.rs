// BAM/CRAM file reader — full implementation using noodles

use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;

use crate::basics::aligned_read::AlignedRead;
use crate::basics::cigar_string::{CigarFlag, CigarOperation, CigarString};
use crate::basics::genomic_region::GenomicRegion;

use super::read_manager::{ReadContainer, ReadReaderImpl, SampleName};

pub struct BamReader {
    path: PathBuf,
    samples: Vec<SampleName>,
    reads_by_sample: HashMap<String, Vec<AlignedRead>>,
}

impl BamReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, String> {
        let path = path.as_ref().to_path_buf();

        // ── First pass: read header ────────────────────────────────────────────
        let file = File::open(&path)
            .map_err(|e| format!("cannot open BAM '{}': {}", path.display(), e))?;
        let mut reader = bam::io::Reader::new(BufReader::new(file));
        let header = reader
            .read_header()
            .map_err(|e| format!("cannot read BAM header '{}': {}", path.display(), e))?;

        // Reference sequence names ordered by BAM ref ID (index = ref_id)
        let ref_names: Vec<String> = header
            .reference_sequences()
            .keys()
            .map(|k| k.to_string())
            .collect();

        // Use filename stem as the default sample name.
        // (Proper @RG SM: parsing is version-dependent; filename is always available.)
        let default_sample = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("SAMPLE")
            .to_string();
        let samples = vec![default_sample.clone()];

        // ── Second pass: read all records into memory ──────────────────────────
        let file2 = File::open(&path)
            .map_err(|e| format!("cannot open BAM '{}': {}", path.display(), e))?;
        let mut reader2 = bam::io::Reader::new(BufReader::new(file2));
        reader2.read_header().map_err(|e| e.to_string())?;

        let mut reads_by_sample: HashMap<String, Vec<AlignedRead>> = HashMap::new();
        let mut n_skipped = 0usize;

        for result in reader2.records() {
            let record = result.map_err(|e| format!("BAM record error: {}", e))?;

            // ── Filter: skip unmapped, duplicates, QC-fail, secondary, supplementary
            let flags = record.flags();
            if flags.is_unmapped()
                || flags.is_duplicate()
                || flags.is_qc_fail()
                || flags.is_secondary()
                || flags.is_supplementary()
            {
                n_skipped += 1;
                continue;
            }

            // ── Contig name (lazy-decoded: reference_sequence_id returns Option<Result<usize>>)
            let ref_id = match record.reference_sequence_id() {
                Some(Ok(id)) => id,
                Some(Err(_)) | None => { n_skipped += 1; continue; }
            };
            let contig = match ref_names.get(ref_id) {
                Some(name) => name.clone(),
                None => { n_skipped += 1; continue; }
            };

            // ── Alignment start (lazy-decoded: alignment_start returns Option<Result<Position>>)
            // noodles Position is 1-based; we use 0-based internally.
            let aln_start: u32 = match record.alignment_start() {
                Some(Ok(pos)) => usize::from(pos) as u32 - 1,
                Some(Err(_)) | None => { n_skipped += 1; continue; }
            };

            // ── CIGAR — compute reference length on the fly
            let mut cigar_ops: CigarString = Vec::new();
            let mut ref_consumed: u32 = 0;
            for op_res in record.cigar().iter() {
                let op = match op_res {
                    Ok(o) => o,
                    Err(_) => { n_skipped += 1; break; }
                };
                let flag = match op.kind() {
                    CigarKind::Match            => CigarFlag::AlignmentMatch,
                    CigarKind::Insertion        => CigarFlag::Insertion,
                    CigarKind::Deletion         => CigarFlag::Deletion,
                    CigarKind::Skip             => CigarFlag::Skipped,
                    CigarKind::SoftClip         => CigarFlag::SoftClipped,
                    CigarKind::HardClip         => CigarFlag::HardClipped,
                    CigarKind::Pad              => CigarFlag::Padding,
                    CigarKind::SequenceMatch    => CigarFlag::SequenceMatch,
                    CigarKind::SequenceMismatch => CigarFlag::Substitution,
                };
                let len = op.len() as u32;
                match flag {
                    CigarFlag::AlignmentMatch
                    | CigarFlag::SequenceMatch
                    | CigarFlag::Substitution
                    | CigarFlag::Deletion
                    | CigarFlag::Skipped => ref_consumed += len,
                    _ => {}
                }
                cigar_ops.push(CigarOperation::new(len, flag));
            }

            if cigar_ops.is_empty() { n_skipped += 1; continue; }

            let aln_end = aln_start + ref_consumed;
            if aln_end <= aln_start { n_skipped += 1; continue; }

            let region = match GenomicRegion::new(&contig, aln_start, aln_end) {
                Ok(r) => r,
                Err(_) => { n_skipped += 1; continue; }
            };

            // ── Sequence: noodles Base → ASCII byte
            let sequence: Vec<u8> = record.sequence().iter().map(|b| u8::from(b)).collect();
            if sequence.is_empty() { n_skipped += 1; continue; }

            // ── Quality scores: QualityScores doesn't expose iter() in 0.65;
            //    use AsRef<[u8]> to access the raw byte slice.
            //    Bind to a local first so the temporary outlives the borrow.
            let qs = record.quality_scores();
            let raw_quals: &[u8] = qs.as_ref();
            let qualities: Vec<u8> = if raw_quals.iter().all(|&q| q == 0xFF) {
                // Missing quality scores encoded as 0xFF; use default Q30
                vec![30u8; sequence.len()]
            } else {
                let mut q = raw_quals.to_vec();
                if q.len() < sequence.len() {
                    q.resize(sequence.len(), 30);
                }
                q
            };

            // ── Read name
            let name: String = record
                .name()
                .map(|n| String::from_utf8_lossy(n.as_ref()).to_string())
                .unwrap_or_else(|| "unknown".to_string());

            // ── Mapping quality
            let mapq: u8 = record.mapping_quality().map(u8::from).unwrap_or(0);

            // ── Flags as u16
            let flags_u16 = u16::from(record.flags());

            let aligned_read = AlignedRead::new(
                name, region, sequence, qualities, cigar_ops, mapq, flags_u16,
            );

            reads_by_sample
                .entry(default_sample.clone())
                .or_default()
                .push(aligned_read);
        }

        let total: usize = reads_by_sample.values().map(|v| v.len()).sum();
        eprintln!(
            "[bam] sample='{}' loaded {} reads ({} skipped) from '{}'",
            default_sample, total, n_skipped, path.display()
        );

        Ok(BamReader { path, samples, reads_by_sample })
    }

    pub fn path(&self) -> &Path { &self.path }

    pub fn total_reads(&self) -> usize {
        self.reads_by_sample.values().map(|v| v.len()).sum()
    }
}

impl ReadReaderImpl for BamReader {
    fn samples(&self) -> Vec<SampleName> {
        self.samples.clone()
    }

    fn fetch(&self, sample: &str, region: &GenomicRegion) -> ReadContainer {
        let reads = match self.reads_by_sample.get(sample) {
            Some(v) => v,
            None => return Vec::new(),
        };
        reads
            .iter()
            .filter(|r| {
                r.mapped_region().contig_name() == region.contig_name()
                    && r.mapped_region().begin() < region.end()
                    && r.mapped_region().end() > region.begin()
            })
            .cloned()
            .collect()
    }
}
