// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Converted from C++ to Rust

use std::collections::HashMap;
use crate::io::variant::vcf_record::VcfRecord;
use crate::basics::genomic_region::GenomicRegion;
use crate::core::types::variant::Variant;
use crate::core::types::haplotype::Haplotype;
use crate::core::types::genotype::{Genotype, make_genotypes};
use crate::core::models::haplotype_likelihood::{HaplotypeLikelihoodArray, BaseQualityLikelihoodModel};
use crate::core::models::genotype::individual_model::{IndividualModel, GenotypeInferenceResult};
use crate::core::models::genotype::genotype_prior_model::HardyWeinbergGenotypePriorModel;
use crate::core::tools::vargen::cigar_scanner::CigarScanner;
use crate::core::tools::vargen::variant_generator::VariantGenerator;
use super::caller::{Caller, CallerEnvironment, CallerOptions, CallBuffer};

pub struct IndividualCaller {
    options: CallerOptions,
}

impl IndividualCaller {
    pub fn new(options: CallerOptions) -> Self {
        IndividualCaller { options }
    }
}

impl Caller for IndividualCaller {
    fn name(&self) -> &str { "individual" }

    fn call_variants(&self, env: &CallerEnvironment) -> Result<CallBuffer, String> {
        let mut calls = CallBuffer::new();

        // ── 1. Collect all reads for candidate generation ─────────────────────
        let all_reads: Vec<_> = env.reads.values()
            .flat_map(|r| r.iter().cloned())
            .collect();

        if all_reads.is_empty() { return Ok(calls); }

        // ── 2. Candidate variant generation via CIGAR-pileup scan ─────────────
        let scanner    = CigarScanner::new(self.options.min_base_quality);
        let candidates = scanner.generate(&all_reads, &env.region, env.reference);

        if candidates.is_empty() { return Ok(calls); }

        // ── 3. Fetch reference sequence for the calling region ────────────────
        let ref_seq = env.reference.fetch_sequence(&env.region)
            .map_err(|e| e.to_string())?;

        // ── 4. Build haplotypes: reference + one per unique candidate ─────────
        let ref_haplotype = Haplotype::new(env.region.clone(), ref_seq.clone());

        let mut alt_pairs: Vec<(Variant, Haplotype)> = candidates.iter()
            .map(|v| {
                let h = build_alt_haplotype(&env.region, &ref_seq, v);
                (v.clone(), h)
            })
            .collect();

        // Drop alts identical to ref (degenerate indels etc.).
        alt_pairs.retain(|(_, h)| h != &ref_haplotype);

        // Deduplicate by haplotype hash (keeps first occurrence per hash).
        let mut seen_hashes = std::collections::HashSet::new();
        alt_pairs.retain(|(_, h)| seen_hashes.insert(h.get_hash()));

        // Respect max_haplotypes cap (one slot reserved for reference).
        alt_pairs.truncate(self.options.max_haplotypes.saturating_sub(1));

        if alt_pairs.is_empty() { return Ok(calls); }

        let haplotypes: Vec<Haplotype> = std::iter::once(ref_haplotype.clone())
            .chain(alt_pairs.iter().map(|(_, h)| h.clone()))
            .collect();

        // ── 5. Per-sample Bayesian genotype calling ───────────────────────────
        for (sample, reads) in &env.reads {
            if reads.len() < self.options.min_read_depth { continue; }

            // P(read | haplotype) using per-base quality scores.
            let model = BaseQualityLikelihoodModel::new();
            let mut likelihood_array = HaplotypeLikelihoodArray::new();
            likelihood_array.populate(&haplotypes, reads, &model);

            // Enumerate diploid genotypes (combinations-with-repetition).
            let genotypes = make_genotypes(&haplotypes, self.options.ploidy);
            if genotypes.is_empty() { continue; }

            // Hardy-Weinberg prior: P(G) = C(n,k)·θ^k·(1−θ)^(n−k).
            let prior_model = HardyWeinbergGenotypePriorModel::new(
                self.options.snp_heterozygosity,
                self.options.indel_heterozygosity,
                ref_haplotype.get_hash(),
            );

            let inference = IndividualModel::new(&prior_model);
            let result    = inference.evaluate(&genotypes, &likelihood_array);

            // ── Variant-level quality filter ──────────────────────────────────
            let var_post   = variant_posterior(&result, &ref_haplotype);
            let qual_phred = prob_to_phred(var_post);
            if qual_phred < self.options.min_variant_quality { continue; }

            let map_genotype  = &genotypes[result.map_genotype_index];
            let map_posterior = result.posteriors[result.map_genotype_index].1;

            // Skip homozygous-reference calls.
            if is_hom_ref(map_genotype, &ref_haplotype) { continue; }

            // ── 6. One VcfRecord per unique non-ref haplotype in MAP genotype ──
            let map_unique = map_genotype.unique_haplotypes();

            for &alt_hap in &map_unique {
                if alt_hap == &ref_haplotype { continue; }

                // Find originating variant for this alt haplotype.
                let variant = match alt_pairs.iter()
                    .find(|(_, h)| h == alt_hap)
                    .map(|(v, _)| v)
                {
                    Some(v) => v,
                    None    => continue,
                };

                // Build properly padded VCF alleles.
                let (vcf_pos, ref_str, alt_str) =
                    make_vcf_alleles(variant, &ref_seq, env.region.begin());

                // GT string, e.g. "0/1" or "1/1".
                let all_alts: Vec<&Haplotype> = alt_pairs.iter().map(|(_, h)| h).collect();
                let gt = gt_string(map_genotype, &ref_haplotype, &all_alts);

                // GQ: Phred(P(wrong MAP genotype)).
                let gq = prob_to_phred(map_posterior) as u32;

                // DP and per-allele depth.
                let dp = reads.len();
                let (ref_depth, alt_depth) = allelic_depths(reads, variant);

                // VcfRecord region (0-based, length = len(REF)).
                let record_region = GenomicRegion::new(
                    variant.mapped_region().contig_name(),
                    vcf_pos,
                    vcf_pos + ref_str.len() as u32,
                ).map_err(|e| e.to_string())?;

                let mut sample_fields: HashMap<String, Vec<String>> = HashMap::new();
                sample_fields.insert("GT".to_string(), vec![gt]);
                sample_fields.insert("GQ".to_string(), vec![gq.to_string()]);
                sample_fields.insert("DP".to_string(), vec![dp.to_string()]);
                sample_fields.insert("AD".to_string(),
                    vec![ref_depth.to_string(), alt_depth.to_string()]);

                let mut sample_map: HashMap<String, HashMap<String, Vec<String>>> = HashMap::new();
                sample_map.insert(sample.clone(), sample_fields);

                let record = VcfRecord::new(
                    record_region,
                    ".".to_string(),
                    ref_str,
                    vec![alt_str],
                    Some(qual_phred as f32),
                    vec!["PASS".to_string()],
                    HashMap::new(),
                ).with_genotypes(
                    vec!["GT".to_string(), "GQ".to_string(),
                         "DP".to_string(), "AD".to_string()],
                    sample_map,
                );

                calls.push(record);
            }
        }

        Ok(calls)
    }
}

// ── Private helpers ────────────────────────────────────────────────────────────

/// Apply `variant` to `ref_seq` (which covers `region`) to produce the alt
/// haplotype sequence, then wrap it in a `Haplotype`.
fn build_alt_haplotype(region: &GenomicRegion, ref_seq: &[u8], variant: &Variant) -> Haplotype {
    let region_begin = region.begin();
    let var_begin    = variant.mapped_region().begin();
    let var_end      = variant.mapped_region().end();

    let off_begin = var_begin.saturating_sub(region_begin) as usize;
    let off_end   = var_end.saturating_sub(region_begin) as usize;
    let off_begin = off_begin.min(ref_seq.len());
    let off_end   = off_end.min(ref_seq.len());

    let mut alt_seq = Vec::with_capacity(ref_seq.len() + variant.alt_sequence_size() + 4);
    alt_seq.extend_from_slice(&ref_seq[..off_begin]);
    alt_seq.extend_from_slice(variant.alt_sequence());
    alt_seq.extend_from_slice(&ref_seq[off_end..]);

    Haplotype::new(region.clone(), alt_seq)
}

/// Convert a posterior probability to a Phred-scaled quality score.
/// QUAL = −10·log₁₀(1 − P(variant)).
fn prob_to_phred(p: f64) -> f64 {
    let p_wrong = (1.0 - p).clamp(1e-30, 1.0);
    -10.0 * p_wrong.log10()
}

/// P(at least one alt allele) = 1 − P(homozygous reference).
fn variant_posterior(result: &GenotypeInferenceResult, ref_haplotype: &Haplotype) -> f64 {
    let hom_ref_prob: f64 = result.posteriors.iter()
        .filter(|(g, _)| g.haplotypes().iter().all(|h| h == ref_haplotype))
        .map(|(_, p)| *p)
        .sum();
    (1.0_f64 - hom_ref_prob).clamp(0.0, 1.0)
}

/// True when the genotype is homozygous for the reference haplotype.
fn is_hom_ref(genotype: &Genotype, ref_haplotype: &Haplotype) -> bool {
    genotype.is_homozygous()
        && genotype.haplotypes().first().map(|h| h == ref_haplotype).unwrap_or(true)
}

/// Build the unphased GT string (e.g. "0/1", "1/1") for the MAP genotype.
/// Allele 0 = reference; allele i+1 = `all_alts[i]`.
fn gt_string(
    map_genotype: &Genotype,
    ref_haplotype: &Haplotype,
    all_alts: &[&Haplotype],
) -> String {
    let mut indices: Vec<usize> = map_genotype.haplotypes().iter()
        .map(|h| {
            if h == ref_haplotype {
                0
            } else {
                all_alts.iter()
                    .position(|ah| *ah == h)
                    .map(|i| i + 1)
                    .unwrap_or(0)
            }
        })
        .collect();
    indices.sort_unstable();
    indices.iter().map(|i| i.to_string()).collect::<Vec<_>>().join("/")
}

/// Build VCF-spec allele strings, adding an anchor base for indels/empty alleles.
/// Returns `(vcf_pos_0based, ref_allele_str, alt_allele_str)`.
fn make_vcf_alleles(
    variant: &Variant,
    ref_seq: &[u8],
    region_begin: u32,
) -> (u32, String, String) {
    let ref_allele = variant.ref_sequence();
    let alt_allele = variant.alt_sequence();
    let var_begin  = variant.mapped_region().begin();

    if ref_allele.is_empty() || alt_allele.is_empty() {
        // Insertion or deletion: VCF requires a non-empty anchor base.
        if var_begin > region_begin {
            let anchor_idx = (var_begin - 1 - region_begin) as usize;
            if let Some(&anchor_byte) = ref_seq.get(anchor_idx) {
                let anchor_ch = char::from(anchor_byte).to_string();
                let ref_str   = anchor_ch.clone() + &String::from_utf8_lossy(ref_allele);
                let alt_str   = anchor_ch          + &String::from_utf8_lossy(alt_allele);
                return (var_begin - 1, ref_str, alt_str);
            }
        }
        // Fallback: variant is at the very start of the region.
        let ref_str = "N".to_string() + &String::from_utf8_lossy(ref_allele);
        let alt_str = "N".to_string() + &String::from_utf8_lossy(alt_allele);
        (var_begin.saturating_sub(1), ref_str, alt_str)
    } else {
        // SNV or MNV: no anchor needed.
        (
            var_begin,
            String::from_utf8_lossy(ref_allele).to_string(),
            String::from_utf8_lossy(alt_allele).to_string(),
        )
    }
}

/// Count reads supporting the reference and alt alleles at the variant position.
/// Best suited for SNVs; approximate for indels.
fn allelic_depths(
    reads: &[crate::basics::aligned_read::AlignedRead],
    variant: &Variant,
) -> (usize, usize) {
    let pos      = variant.mapped_region().begin();
    let ref_base = variant.ref_sequence().first().copied();
    let alt_base = variant.alt_sequence().first().copied();

    let mut ref_depth = 0usize;
    let mut alt_depth = 0usize;

    for read in reads {
        let r = read.mapped_region();
        if pos < r.begin() || pos >= r.end() { continue; }
        let idx = (pos - r.begin()) as usize;
        if let Some(&base) = read.sequence().get(idx) {
            if Some(base) == alt_base {
                alt_depth += 1;
            } else if Some(base) == ref_base {
                ref_depth += 1;
            }
        }
    }

    (ref_depth, alt_depth)
}

// ── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use crate::basics::aligned_read::AlignedRead;
    use crate::basics::cigar_string::parse_cigar;
    use crate::basics::genomic_region::GenomicRegion;
    use crate::core::callers::caller::{CallerEnvironment, CallerOptions};
    use crate::io::reference::reference_genome::{ReferenceGenome, ReferenceReader, GeneticSequence};
    use crate::basics::contig_region::Size;

    // ── Minimal in-memory reference ───────────────────────────────────────────

    struct MockReference {
        sequence: Vec<u8>,
        contig:   String,
    }

    impl ReferenceReader for MockReference {
        fn name(&self) -> &str { "mock" }
        fn contig_names(&self) -> Vec<String> { vec![self.contig.clone()] }
        fn contig_size(&self, _: &str) -> Option<Size> {
            Some(self.sequence.len() as Size)
        }
        fn fetch_sequence(&self, region: &GenomicRegion) -> Result<GeneticSequence, String> {
            let begin = region.begin() as usize;
            let end   = region.end()   as usize;
            if end > self.sequence.len() {
                return Err(format!(
                    "region {}..{} out of bounds (len={})", begin, end, self.sequence.len()
                ));
            }
            Ok(self.sequence[begin..end].to_vec())
        }
    }

    fn make_ref(contig: &str, seq: &[u8]) -> ReferenceGenome {
        ReferenceGenome::new(Box::new(MockReference {
            sequence: seq.to_vec(),
            contig: contig.to_string(),
        }))
    }

    fn make_read(contig: &str, begin: u32, seq: &[u8], qual: u8, cigar: &str) -> AlignedRead {
        let end    = begin + seq.len() as u32;
        let region = GenomicRegion::new(contig, begin, end).unwrap();
        let quals  = vec![qual; seq.len()];
        AlignedRead::new(
            "read".to_string(),
            region,
            seq.to_vec(),
            quals,
            parse_cigar(cigar).unwrap(),
            60,
            0x1,
        )
    }

    // ── build_alt_haplotype ───────────────────────────────────────────────────

    #[test]
    fn build_alt_haplotype_snv() {
        let region  = GenomicRegion::new("chr1", 0, 8).unwrap();
        let ref_seq = b"ACGTACGT";
        // T→G at position 3: ref[0..3]=ACG + G + ref[4..8]=ACGT → ACGGACGT (8 bases)
        let variant = Variant::from_parts("chr1", 3, b"T".to_vec(), b"G".to_vec()).unwrap();
        let alt     = build_alt_haplotype(&region, ref_seq, &variant);
        assert_eq!(alt.sequence(), b"ACGGACGT".as_ref());
    }

    #[test]
    fn build_alt_haplotype_insertion() {
        let region  = GenomicRegion::new("chr1", 0, 8).unwrap();
        let ref_seq = b"ACGTACGT";
        // Insert "TT" before position 4 (ref allele empty)
        let variant = Variant::from_parts("chr1", 4, b"".to_vec(), b"TT".to_vec()).unwrap();
        let alt     = build_alt_haplotype(&region, ref_seq, &variant);
        // ref[0..4]=ACGT + TT + ref[4..]=ACGT → ACGTTTACGT
        assert_eq!(alt.sequence(), b"ACGTTTACGT".as_ref());
    }

    #[test]
    fn build_alt_haplotype_deletion() {
        let region  = GenomicRegion::new("chr1", 0, 8).unwrap();
        let ref_seq = b"ACGTACGT";
        // Delete GT at positions 2-3 (ref = GT, alt empty)
        let variant = Variant::from_parts("chr1", 2, b"GT".to_vec(), b"".to_vec()).unwrap();
        let alt     = build_alt_haplotype(&region, ref_seq, &variant);
        // ref[0..2]=AC + "" + ref[4..]=ACGT → ACACGT
        assert_eq!(alt.sequence(), b"ACACGT".as_ref());
    }

    // ── prob_to_phred ─────────────────────────────────────────────────────────

    #[test]
    fn prob_to_phred_near_one_gives_high_score() {
        assert!(prob_to_phred(0.999) > 25.0);
    }

    #[test]
    fn prob_to_phred_half_gives_three() {
        assert!((prob_to_phred(0.5) - 3.0103).abs() < 0.01);
    }

    // ── gt_string ─────────────────────────────────────────────────────────────

    #[test]
    fn gt_string_hom_ref() {
        let region = GenomicRegion::new("chr1", 0, 4).unwrap();
        let ref_h  = Haplotype::new(region.clone(), b"ACGT".to_vec());
        let gt     = Genotype::new(vec![ref_h.clone(), ref_h.clone()]);
        let alts: Vec<&Haplotype> = vec![];
        assert_eq!(gt_string(&gt, &ref_h, &alts), "0/0");
    }

    #[test]
    fn gt_string_het() {
        let region = GenomicRegion::new("chr1", 0, 4).unwrap();
        let ref_h  = Haplotype::new(region.clone(), b"ACGT".to_vec());
        let alt_h  = Haplotype::new(region.clone(), b"ACTT".to_vec());
        let gt     = Genotype::new(vec![ref_h.clone(), alt_h.clone()]);
        let alts   = vec![&alt_h];
        assert_eq!(gt_string(&gt, &ref_h, &alts), "0/1");
    }

    #[test]
    fn gt_string_hom_alt() {
        let region = GenomicRegion::new("chr1", 0, 4).unwrap();
        let ref_h  = Haplotype::new(region.clone(), b"ACGT".to_vec());
        let alt_h  = Haplotype::new(region.clone(), b"ACTT".to_vec());
        let gt     = Genotype::new(vec![alt_h.clone(), alt_h.clone()]);
        let alts   = vec![&alt_h];
        assert_eq!(gt_string(&gt, &ref_h, &alts), "1/1");
    }

    // ── make_vcf_alleles ──────────────────────────────────────────────────────

    #[test]
    fn vcf_alleles_snv_no_anchor() {
        let ref_seq = b"ACGTACGT";
        let variant = Variant::from_parts("chr1", 3, b"T".to_vec(), b"G".to_vec()).unwrap();
        let (pos, r, a) = make_vcf_alleles(&variant, ref_seq, 0);
        assert_eq!(pos, 3);
        assert_eq!(r, "T");
        assert_eq!(a, "G");
    }

    #[test]
    fn vcf_alleles_insertion_uses_anchor() {
        // anchor base at pos 3 = 'T'
        let ref_seq = b"ACGTACGT";
        let variant = Variant::from_parts("chr1", 4, b"".to_vec(), b"TT".to_vec()).unwrap();
        let (pos, r, a) = make_vcf_alleles(&variant, ref_seq, 0);
        assert_eq!(pos, 3);
        assert_eq!(r, "T");
        assert_eq!(a, "TTT");
    }

    #[test]
    fn vcf_alleles_deletion_uses_anchor() {
        // anchor base at pos 1 = 'C'
        let ref_seq = b"ACGTACGT";
        let variant = Variant::from_parts("chr1", 2, b"GT".to_vec(), b"".to_vec()).unwrap();
        let (pos, r, a) = make_vcf_alleles(&variant, ref_seq, 0);
        assert_eq!(pos, 1);
        assert_eq!(r, "CGT");
        assert_eq!(a, "C");
    }

    // ── End-to-end calling ────────────────────────────────────────────────────

    /// Reference: chr1 0-20  = ACGTACGTACGTACGTACGT (20 bases).
    ///   Position 9 (0-based) = C.
    /// 8 reads carry T at position 9 (alt), 2 carry C (ref).
    /// Expected: at least one call at VCF POS=10 with REF=C ALT=T.
    #[test]
    fn calls_snv_het_from_reads() {
        let ref_bases: &[u8] = b"ACGTACGTACGTACGTACGT"; // pos 9 = C
        let reference = make_ref("chr1", ref_bases);
        let region    = GenomicRegion::new("chr1", 0, 20).unwrap();

        //                        0123456789...
        let alt_seq:  &[u8] = b"ACGTACGTATGTACGTACGT"; // pos 9 = T

        let mut reads = Vec::new();
        for _ in 0..8 {
            reads.push(make_read("chr1", 0, alt_seq,   30, "20M"));
        }
        for _ in 0..2 {
            reads.push(make_read("chr1", 0, ref_bases, 30, "20M"));
        }

        let mut sample_reads = HashMap::new();
        sample_reads.insert("SAMPLE".to_string(), reads);

        let env = CallerEnvironment {
            reference: &reference,
            reads: sample_reads,
            region,
        };

        let opts   = CallerOptions { min_variant_quality: 1.0, ..CallerOptions::default() };
        let caller = IndividualCaller::new(opts);
        let calls  = caller.call_variants(&env).unwrap();

        assert!(!calls.is_empty(), "expected at least one variant call");

        let record = &calls[0];
        assert_eq!(record.chrom(), "chr1");
        assert_eq!(record.pos(), 10); // 1-based POS for 0-based position 9
        assert_eq!(record.ref_allele(), "C");
        assert_eq!(record.alt(), &["T"]);

        let gt_vals = record.get_sample_value("SAMPLE", "GT").unwrap();
        let gt_val  = &gt_vals[0];
        assert!(
            gt_val == "0/1" || gt_val == "1/1",
            "unexpected GT: {}", gt_val
        );
    }

    /// All reads match reference → no variant call.
    #[test]
    fn no_call_when_reads_match_reference() {
        let ref_bases = b"ACGTACGTACGT";
        let reference = make_ref("chr1", ref_bases);
        let region    = GenomicRegion::new("chr1", 0, 12).unwrap();

        let reads: Vec<AlignedRead> = (0..5)
            .map(|_| make_read("chr1", 0, ref_bases, 30, "12M"))
            .collect();

        let mut sample_reads = HashMap::new();
        sample_reads.insert("SAMPLE".to_string(), reads);

        let env = CallerEnvironment { reference: &reference, reads: sample_reads, region };

        let calls = IndividualCaller::new(CallerOptions::default())
            .call_variants(&env).unwrap();

        assert!(calls.is_empty(), "expected no calls for a perfect-match pileup");
    }

    /// Empty read set → no call.
    #[test]
    fn no_call_on_empty_reads() {
        let ref_bases = b"ACGTACGT";
        let reference = make_ref("chr1", ref_bases);
        let region    = GenomicRegion::new("chr1", 0, 8).unwrap();

        let mut sample_reads = HashMap::new();
        sample_reads.insert("SAMPLE".to_string(), Vec::<AlignedRead>::new());

        let env = CallerEnvironment { reference: &reference, reads: sample_reads, region };

        let calls = IndividualCaller::new(CallerOptions::default())
            .call_variants(&env).unwrap();

        assert!(calls.is_empty());
    }
}
