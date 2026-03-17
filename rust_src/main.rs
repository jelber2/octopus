// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
// Rust port of Octopus variant caller

mod basics;
mod concepts;
mod utils;
mod exceptions;
mod config;
mod logging;
mod io;
mod core;
mod containers;
mod readpipe;

use std::path::PathBuf;
use std::process;
use clap::{Parser, ArgGroup};

use crate::io::reference::fasta::FastaReader;
use crate::io::reference::reference_genome::{ReferenceGenome, get_all_contig_regions};
use crate::io::read::read_manager::ReadManager;
use crate::io::variant::vcf_header::VcfHeader;
use crate::io::variant::vcf_writer::VcfWriter;
use crate::core::callers::caller::{Caller, CallerEnvironment, CallerOptions};
use crate::core::callers::individual_caller::IndividualCaller;
use crate::basics::genomic_region::GenomicRegion;

const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Octopus — a mapping-based haplotype-aware variant caller.
///
/// Octopus calls SNVs, indels, and complex variants across 6 calling models:
/// individual, population, trio, cancer, polyclone, and cell.
#[derive(Parser, Debug)]
#[command(
    name = "octopus",
    version = VERSION,
    about = "A mapping-based haplotype-aware variant caller",
    long_about = None,
    group(
        ArgGroup::new("input_reads")
            .required(true)
            .args(["reads", "reads_file"])
    ),
)]
struct OctopusArgs {
    // ── General ──────────────────────────────────────────────────────────────

    /// Config file to populate command line options
    #[arg(long, value_name = "FILE")]
    config: Option<PathBuf>,

    /// Create log file for debugging [default filename: octopus_debug.log]
    #[arg(long, value_name = "FILE", num_args = 0..=1,
          default_missing_value = "octopus_debug.log")]
    debug: Option<PathBuf>,

    /// Create very verbose log file for debugging [default filename: octopus_trace.log]
    #[arg(long, value_name = "FILE", num_args = 0..=1,
          default_missing_value = "octopus_trace.log")]
    trace: Option<PathBuf>,

    /// Set the working directory
    #[arg(short = 'w', long, value_name = "DIR")]
    working_directory: Option<PathBuf>,

    /// Maximum number of threads to use (0 = unlimited)
    #[arg(long, visible_alias = "tentacles", short = 'j',
          value_name = "N", num_args = 0..=1, default_missing_value = "0")]
    threads: Option<usize>,

    /// Maximum memory for cached reference sequence (e.g. 500MB)
    #[arg(short = 'X', long, value_name = "MEM", default_value = "500MB")]
    max_reference_cache_memory: String,

    /// Target memory limit for buffered read data (e.g. 6GB)
    #[arg(short = 'B', long, value_name = "MEM", default_value = "6GB")]
    target_read_buffer_memory: String,

    /// Maximum number of BAM/CRAM files open simultaneously
    #[arg(long, value_name = "N", default_value = "250")]
    max_open_read_files: usize,

    /// Indexed FASTA reference genome [REQUIRED]
    #[arg(short = 'R', long, value_name = "FILE", required = true)]
    reference: PathBuf,

    /// Indexed BAM/CRAM read files to analyse
    #[arg(short = 'I', long, value_name = "FILE", num_args = 1..)]
    reads: Vec<PathBuf>,

    /// File(s) containing lists of BAM/CRAM paths (one per line), or BAM/CRAM files directly
    #[arg(short = 'i', long, value_name = "FILE", num_args = 1..)]
    reads_file: Vec<PathBuf>,

    /// Regions to analyse, e.g. chr1 or chr1:100-200
    #[arg(short = 'T', long, value_name = "REGION", num_args = 1..)]
    regions: Vec<String>,

    /// File with regions to analyse (one per line)
    #[arg(short = 't', long, value_name = "FILE")]
    regions_file: Option<PathBuf>,

    /// Regions to skip
    #[arg(short = 'K', long, value_name = "REGION", num_args = 1..)]
    skip_regions: Vec<String>,

    /// File with regions to skip (one per line)
    #[arg(short = 'k', long, value_name = "FILE")]
    skip_regions_file: Option<PathBuf>,

    /// Assume 1-based (rather than 0-based) coordinates for --regions / --skip-regions
    #[arg(long, default_value = "false")]
    one_based_indexing: bool,

    /// Subset of sample names to analyse
    #[arg(short = 'S', long, value_name = "SAMPLE", num_args = 1..)]
    samples: Vec<String>,

    /// File of sample names to analyse (one per line)
    #[arg(short = 's', long, value_name = "FILE")]
    samples_file: Option<PathBuf>,

    /// Ignore contigs not mapped in the read files
    #[arg(long, default_value = "false")]
    ignore_unmapped_contigs: bool,

    /// PED pedigree file
    #[arg(long, value_name = "FILE")]
    pedigree: Option<PathBuf>,

    /// Output file for variant calls (stdout if omitted)
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Only report call sites — drop sample genotype information
    #[arg(long, default_value = "false")]
    sites_only: bool,

    /// Output realigned evidence BAM
    #[arg(long, value_name = "FILE")]
    bamout: Option<PathBuf>,

    /// Type of realigned BAM to output [mini, full]
    #[arg(long, value_name = "TYPE", default_value = "mini")]
    bamout_type: String,

    /// Improve runtime at some cost to accuracy (disables several features)
    #[arg(long, default_value = "false")]
    fast: bool,

    /// Like --fast but faster still
    #[arg(long, default_value = "false")]
    very_fast: bool,

    // ── Read preprocessing ────────────────────────────────────────────────

    /// Disable all read preprocessing
    #[arg(long, default_value = "false")]
    disable_read_preprocessing: bool,

    /// Cap all base qualities to this value
    #[arg(long, value_name = "Q")]
    max_base_quality: Option<u8>,

    /// Mask read tail bases with base quality below this
    #[arg(long, value_name = "Q")]
    mask_low_quality_tails: Option<u8>,

    /// Disable adapter detection and masking
    #[arg(long, default_value = "false")]
    disable_adapter_masking: bool,

    /// Disable read-segment overlap masking
    #[arg(long, default_value = "false")]
    disable_overlap_masking: bool,

    /// Minimum read mapping quality to use a read for calling
    #[arg(long, value_name = "Q", default_value = "5")]
    min_mapping_quality: u8,

    /// Base quality threshold for --min-good-bases / --min-good-base-fraction
    #[arg(long, value_name = "Q", default_value = "20")]
    good_base_quality: u8,

    /// Minimum number of bases above --good-base-quality before a read is used
    #[arg(long, value_name = "N", default_value = "20")]
    min_good_bases: u32,

    /// Pass QC-failed reads through
    #[arg(long, default_value = "false")]
    allow_qc_fails: bool,

    /// Filter reads shorter than this
    #[arg(long, value_name = "N")]
    min_read_length: Option<u32>,

    /// Filter reads longer than this
    #[arg(long, value_name = "N", default_value = "10000")]
    max_read_length: u32,

    /// Allow reads flagged as duplicates
    #[arg(long, default_value = "false")]
    allow_marked_duplicates: bool,

    /// Allow secondary alignments
    #[arg(long, default_value = "false")]
    allow_secondary_alignments: bool,

    /// Allow supplementary alignments
    #[arg(long, default_value = "false")]
    allow_supplementary_alignments: bool,

    /// Disable downsampling
    #[arg(long, default_value = "false")]
    disable_downsampling: bool,

    /// Downsample regions with coverage above this
    #[arg(long, value_name = "N", default_value = "1000")]
    downsample_above: u32,

    /// Target coverage after downsampling
    #[arg(long, value_name = "N", default_value = "500")]
    downsample_target: u32,

    // ── Variant discovery ─────────────────────────────────────────────────

    /// Candidate variant discovery protocol [illumina, pacbio]
    #[arg(long, value_name = "PROTOCOL", default_value = "illumina")]
    variant_discovery_mode: String,

    /// Disable all candidate discovery from reads
    #[arg(long, default_value = "false")]
    disable_denovo_variant_discovery: bool,

    /// Disable CIGAR-based pileup candidate generator
    #[arg(long, default_value = "false")]
    disable_pileup_candidate_generator: bool,

    /// Disable local re-assembly candidate generator
    #[arg(long, default_value = "false")]
    disable_assembly_candidate_generator: bool,

    /// VCF files of known variants to include as candidates
    #[arg(short = 'c', long, value_name = "FILE", num_args = 1..)]
    source_candidates: Vec<PathBuf>,

    /// Minimum base quality for pileup candidate generation
    #[arg(long, value_name = "Q", default_value = "20")]
    min_pileup_base_quality: u8,

    /// Minimum reads supporting a variant for it to be a candidate
    #[arg(long, value_name = "N")]
    min_supporting_reads: Option<u32>,

    /// Maximum candidate variant size (reference bases)
    #[arg(long, value_name = "N", default_value = "2000")]
    max_variant_size: u32,

    /// k-mer sizes for local assembly
    #[arg(long, value_name = "K", num_args = 1.., default_values = ["10", "15", "20"])]
    kmer_sizes: Vec<u32>,

    // ── Variant calling (general) ─────────────────────────────────────────

    /// Calling model to use [individual, population, trio, cancer, polyclone, cell]
    #[arg(short = 'C', long, value_name = "MODEL", default_value = "population",
          value_parser = validate_caller_name)]
    caller: String,

    /// Default ploidy for contigs without an explicit ploidy
    #[arg(short = 'P', long, value_name = "N", default_value = "2")]
    organism_ploidy: u32,

    /// Contig-specific ploidies, e.g. Y=1 chrY=1 MT=1
    #[arg(short = 'p', long, value_name = "CONTIG=PLOIDY", num_args = 1..)]
    contig_ploidies: Vec<String>,

    /// Minimum variant posterior (Phred) to emit a call
    #[arg(long, value_name = "PHRED", default_value = "0.1")]
    min_variant_posterior: f64,

    /// Emit reference confidence calls for non-variant positions
    #[arg(long, default_value = "false")]
    refcall: bool,

    /// Germline SNP heterozygosity prior
    #[arg(short = 'z', long, value_name = "PROB", default_value = "0.001")]
    snp_heterozygosity: f64,

    /// Germline indel heterozygosity prior
    #[arg(short = 'y', long, value_name = "PROB", default_value = "0.0001")]
    indel_heterozygosity: f64,

    /// Use uniform genotype priors
    #[arg(long, default_value = "false")]
    use_uniform_genotype_priors: bool,

    /// Maximum number of genotypes to evaluate
    #[arg(long, value_name = "N")]
    max_genotypes: Option<u32>,

    /// Policy for model posterior calculation [all, off, special]
    #[arg(long, value_name = "POLICY", default_value = "all")]
    model_posterior: String,

    /// Ignore read mapping quality in haplotype likelihood calculation
    #[arg(long, default_value = "false")]
    dont_model_mapping_quality: bool,

    /// Sequencing error model for haplotype likelihood
    #[arg(long, value_name = "MODEL", default_value = "PCR-free.HiSeq-2500")]
    sequence_error_model: String,

    /// Read linkage information [none, paired, linked]
    #[arg(long, value_name = "LINKAGE", default_value = "paired")]
    read_linkage: String,

    /// Minimum phase score (Phred) to report phased sites
    #[arg(long, value_name = "PHRED", default_value = "5.0")]
    min_phase_score: f64,

    /// Phasing policy [auto, conservative, aggressive]
    #[arg(long, value_name = "POLICY", default_value = "auto")]
    phasing_policy: String,

    // ── Cancer calling model ──────────────────────────────────────────────

    /// Normal sample names (all others are treated as tumour)
    #[arg(short = 'N', long, value_name = "SAMPLE", num_args = 1..)]
    normal_samples: Vec<String>,

    /// Prior probability of a somatic SNV at any given base
    #[arg(long, value_name = "PROB", default_value = "0.0001")]
    somatic_snv_prior: f64,

    /// Prior probability of a somatic indel at any given position
    #[arg(long, value_name = "PROB", default_value = "0.000001")]
    somatic_indel_prior: f64,

    /// Minimum credible somatic allele frequency to report
    #[arg(long, value_name = "FREQ", default_value = "0.005")]
    min_credible_somatic_frequency: f64,

    /// Minimum somatic posterior (Phred) to emit a somatic call
    #[arg(long, value_name = "PHRED", default_value = "0.5")]
    min_somatic_posterior: f64,

    /// Only emit SOMATIC mutations (suppress germline calls)
    #[arg(long, default_value = "false")]
    somatics_only: bool,

    // ── Trio calling model ────────────────────────────────────────────────

    /// Maternal sample name
    #[arg(short = 'M', long, value_name = "SAMPLE")]
    maternal_sample: Option<String>,

    /// Paternal sample name
    #[arg(short = 'F', long, value_name = "SAMPLE")]
    paternal_sample: Option<String>,

    /// Prior probability of a de novo SNV in the offspring
    #[arg(long, value_name = "PROB", default_value = "1.3e-8")]
    denovo_snv_prior: f64,

    /// Prior probability of a de novo indel in the offspring
    #[arg(long, value_name = "PROB", default_value = "1e-9")]
    denovo_indel_prior: f64,

    /// Minimum de novo posterior (Phred) to emit a de novo call
    #[arg(long, value_name = "PHRED", default_value = "3.0")]
    min_denovo_posterior: f64,

    /// Only emit DENOVO mutations (suppress inherited variant calls)
    #[arg(long, default_value = "false")]
    denovos_only: bool,

    // ── Polyclone calling model ───────────────────────────────────────────

    /// Maximum number of distinct clones to consider
    #[arg(long, value_name = "N", default_value = "5")]
    max_clones: u32,

    /// Minimum expected clone frequency
    #[arg(long, value_name = "FREQ", default_value = "0.01")]
    min_clone_frequency: f64,

    /// Prior probability of each clone in the sample
    #[arg(long, value_name = "PROB", default_value = "0.1")]
    clone_prior: f64,

    // ── Cell calling model ────────────────────────────────────────────────

    /// Maximum number of haplotype copy losses in the phylogeny
    #[arg(long, value_name = "N", default_value = "0")]
    max_copy_loss: u32,

    /// Maximum number of haplotype copy gains in the phylogeny
    #[arg(long, value_name = "N", default_value = "0")]
    max_copy_gain: u32,

    /// Allelic dropout concentration parameter
    #[arg(long, value_name = "ALPHA", default_value = "5.0")]
    dropout_concentration: f64,
}

fn validate_caller_name(s: &str) -> Result<String, String> {
    match s {
        "individual" | "population" | "trio" | "cancer" | "polyclone" | "cell" => Ok(s.to_string()),
        other => Err(format!(
            "Unknown caller '{}'. Must be one of: individual, population, trio, cancer, polyclone, cell",
            other
        )),
    }
}

// ── BAM path resolution ───────────────────────────────────────────────────────

/// Resolve the final list of BAM file paths from CLI arguments.
/// `--reads` paths are used directly.  `--reads-file` paths are either:
///   - a BAM/CRAM file (detected by extension), used directly; or
///   - a text file containing one BAM/CRAM path per line.
fn resolve_bam_paths(args: &OctopusArgs) -> Result<Vec<PathBuf>, String> {
    let mut paths: Vec<PathBuf> = args.reads.clone();

    for list_path in &args.reads_file {
        let ext = list_path
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| e.to_ascii_lowercase());

        if matches!(ext.as_deref(), Some("bam") | Some("cram")) {
            paths.push(list_path.clone());
        } else {
            let content = std::fs::read_to_string(list_path).map_err(|e| {
                format!("cannot read reads-file '{}': {}", list_path.display(), e)
            })?;
            for line in content.lines() {
                let line = line.trim();
                if !line.is_empty() && !line.starts_with('#') {
                    paths.push(PathBuf::from(line));
                }
            }
        }
    }

    Ok(paths)
}

// ── Windowed contig processing ────────────────────────────────────────────────

/// Process one contig in non-overlapping windows of `window_size` bp.
/// Returns the number of variant records emitted.
fn process_contig(
    reference: &ReferenceGenome,
    read_manager: &ReadManager,
    caller: &dyn Caller,
    contig: &str,
    contig_size: u32,
    window_size: u32,
    vcf_writer: &mut VcfWriter,
) -> Result<usize, String> {
    let mut total_calls = 0usize;
    let mut pos = 0u32;

    while pos < contig_size {
        let win_begin = pos;
        let win_end = (pos + window_size).min(contig_size);

        let region = GenomicRegion::new(contig, win_begin, win_end)
            .map_err(|e| format!("{}: {}..{}: {}", contig, win_begin, win_end, e))?;

        let reads = read_manager.fetch(&region);
        let has_reads = reads.values().any(|v| !v.is_empty());

        if has_reads {
            let env = CallerEnvironment {
                reference,
                reads,
                region: region.clone(),
            };

            match caller.call_variants(&env) {
                Ok(calls) => {
                    for record in calls {
                        vcf_writer.write_record(&record)?;
                        total_calls += 1;
                    }
                }
                Err(e) => {
                    eprintln!("[octopus] warning: error on {}:{}-{}: {}", contig, win_begin, win_end, e);
                }
            }
        }

        pos = win_end;
    }

    Ok(total_calls)
}

// ── Main run logic ────────────────────────────────────────────────────────────

fn run(args: OctopusArgs) -> i32 {
    if let Err(e) = validate_args(&args) {
        eprintln!("error: {}", e);
        return 1;
    }

    // ── Banner ────────────────────────────────────────────────────────────────
    eprintln!("Octopus variant caller v{}", VERSION);
    eprintln!();
    eprintln!("  Reference:  {}", args.reference.display());
    if !args.reads.is_empty() {
        eprintln!(
            "  Reads:      {}",
            args.reads.iter().map(|p| p.display().to_string()).collect::<Vec<_>>().join(", ")
        );
    }
    if !args.reads_file.is_empty() {
        eprintln!(
            "  Reads file: {}",
            args.reads_file.iter().map(|p| p.display().to_string()).collect::<Vec<_>>().join(", ")
        );
    }
    if !args.regions.is_empty() {
        eprintln!("  Regions:    {}", args.regions.join(" "));
    }
    eprintln!("  Caller:     {}", args.caller);
    eprintln!("  Ploidy:     {}", args.organism_ploidy);
    if let Some(ref out) = args.output {
        eprintln!("  Output:     {}", out.display());
    } else {
        eprintln!("  Output:     stdout");
    }
    if let Some(threads) = args.threads {
        eprintln!("  Threads:    {}", if threads == 0 { "unlimited".to_string() } else { threads.to_string() });
    }
    eprintln!();

    // ── Caller-specific validation ────────────────────────────────────────────
    match args.caller.as_str() {
        "cancer" => {
            if args.normal_samples.is_empty() {
                eprintln!("error: --caller cancer requires at least one --normal-samples");
                return 1;
            }
        }
        "trio" => {
            if args.maternal_sample.is_none() || args.paternal_sample.is_none() {
                eprintln!("error: --caller trio requires both --maternal-sample (-M) and --paternal-sample (-F)");
                return 1;
            }
        }
        _ => {}
    }

    // ── Open reference genome (FASTA) ─────────────────────────────────────────
    eprintln!("[octopus] Loading reference...");
    let fasta_reader = match FastaReader::new(&args.reference) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("error: cannot open reference '{}': {}", args.reference.display(), e);
            return 1;
        }
    };
    let reference = ReferenceGenome::new(Box::new(fasta_reader));
    eprintln!("[octopus] Reference has {} contigs.", reference.num_contigs());

    // ── Resolve and open BAM files ────────────────────────────────────────────
    let bam_paths = match resolve_bam_paths(&args) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("error: {}", e);
            return 1;
        }
    };

    if bam_paths.is_empty() {
        eprintln!("error: no BAM/CRAM files specified");
        return 1;
    }

    eprintln!("[octopus] Opening {} BAM file(s)...", bam_paths.len());
    let read_manager = match ReadManager::new(bam_paths, args.max_open_read_files) {
        Ok(rm) => rm,
        Err(e) => {
            eprintln!("error: {}", e);
            return 1;
        }
    };
    eprintln!("[octopus] Samples: {:?}", read_manager.samples());

    // ── Build VCF header ──────────────────────────────────────────────────────
    let mut header = VcfHeader::new("VCFv4.3");
    for contig in reference.contig_names() {
        let len = reference.contig_size(contig).map(|s| s as u64);
        header.add_contig(contig, len);
    }
    header.add_filter("PASS", "All filters passed");
    header.add_format("GT", "1", "String", "Genotype");
    header.add_format("GQ", "1", "Integer", "Genotype quality");
    header.add_format("DP", "1", "Integer", "Total read depth at this position");
    header.add_format("AD", "R", "Integer", "Allelic depth (ref, alt)");
    for sample in read_manager.samples() {
        header.add_sample(sample);
    }

    // ── Open output ───────────────────────────────────────────────────────────
    let mut vcf_writer = match &args.output {
        Some(path) => match VcfWriter::to_file(path) {
            Ok(w) => w,
            Err(e) => {
                eprintln!("error: cannot create output '{}': {}", path.display(), e);
                return 1;
            }
        },
        None => VcfWriter::to_stdout(),
    };

    if let Err(e) = vcf_writer.write_header(&header) {
        eprintln!("error: writing VCF header: {}", e);
        return 1;
    }

    // ── Build caller ──────────────────────────────────────────────────────────
    let options = CallerOptions {
        min_variant_quality: args.min_variant_posterior,
        ploidy: args.organism_ploidy as usize,
        snp_heterozygosity: args.snp_heterozygosity,
        indel_heterozygosity: args.indel_heterozygosity,
        min_base_quality: args.min_pileup_base_quality,
        ..CallerOptions::default()
    };

    let caller: Box<dyn Caller> = match args.caller.as_str() {
        "individual" => Box::new(IndividualCaller::new(options)),
        other => {
            eprintln!("error: caller '{}' not yet implemented", other);
            return 1;
        }
    };

    // ── Process regions ───────────────────────────────────────────────────────
    // Window size: 50 kb balances memory and throughput.
    const WINDOW_SIZE: u32 = 50_000;

    let regions_to_process: Vec<GenomicRegion> = get_all_contig_regions(&reference);

    let mut total_calls = 0usize;
    let num_contigs = regions_to_process.len();

    for (idx, region) in regions_to_process.iter().enumerate() {
        let contig = region.contig_name();
        let contig_size = region.end();
        eprintln!(
            "[octopus] Processing contig {}/{}: {} ({} bp)",
            idx + 1, num_contigs, contig, contig_size
        );

        match process_contig(
            &reference,
            &read_manager,
            caller.as_ref(),
            contig,
            contig_size,
            WINDOW_SIZE,
            &mut vcf_writer,
        ) {
            Ok(n) => {
                total_calls += n;
                if n > 0 {
                    eprintln!("[octopus]   {} variant(s) called on {}", n, contig);
                }
            }
            Err(e) => {
                eprintln!("[octopus] error on contig {}: {}", contig, e);
            }
        }
    }

    if let Err(e) = vcf_writer.flush() {
        eprintln!("error: flushing VCF output: {}", e);
        return 1;
    }

    eprintln!();
    eprintln!("[octopus] Done. Total variant calls: {}", total_calls);
    0
}

fn validate_args(args: &OctopusArgs) -> Result<(), String> {
    if !args.reference.exists() {
        return Err(format!("reference file not found: {}", args.reference.display()));
    }
    for bam in &args.reads {
        if !bam.exists() {
            return Err(format!("read file not found: {}", bam.display()));
        }
    }
    for f in &args.reads_file {
        if !f.exists() {
            return Err(format!("reads-file not found: {}", f.display()));
        }
    }
    if let Some(ref regions_file) = args.regions_file {
        if !regions_file.exists() {
            return Err(format!("regions-file not found: {}", regions_file.display()));
        }
    }
    if args.fast && args.very_fast {
        return Err("--fast and --very-fast are mutually exclusive".to_string());
    }
    if args.downsample_target >= args.downsample_above && !args.disable_downsampling {
        return Err(format!(
            "--downsample-target ({}) must be less than --downsample-above ({})",
            args.downsample_target, args.downsample_above
        ));
    }
    Ok(())
}

fn main() {
    env_logger::init();

    let args = OctopusArgs::parse();
    let exit_code = run(args);
    process::exit(exit_code);
}
