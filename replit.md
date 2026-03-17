# Octopus — Rust Conversion

## Overview
Full conversion of the Octopus bioinformatics variant caller from C++ (~154,000 lines,
636 files) to Rust.  Supports 6 calling models: individual, population, trio, cancer,
polyclone, cell.

## Build
```
bash scripts/build_rust.sh   # restores missing sources, cargo build --release, cargo test
```
Or directly (Rust stable must be on PATH):
```
export PATH="$HOME/.nix-profile/bin:$PATH"
cargo build --release
cargo test
```
Current status: **209 tests passing, 0 failing**.

## Project Structure
```
rust_src/
  main.rs                      CLI entry point (all 6 callers, --help)
  basics/
    genomic_region.rs           Genomic coordinate type
    contig_region.rs            Single-contig coordinate
    aligned_read.rs             BAM read representation
    cigar_string.rs             CIGAR parsing and operations
    phred.rs                    Phred quality score utilities
  core/
    types/
      allele.rs                 Allele and NucleotideSequence types
      haplotype.rs              Haplotype (sequence + region + explicit alleles)
      genotype.rs               Genotype + make_genotypes (combinations-with-repetition)
      variant.rs                Variant (SNV/indel), helper predicates
    callers/
      caller.rs                 Caller trait, CallerEnvironment, CallerOptions
      individual_caller.rs  ★  Full Bayesian SNV/indel calling pipeline
      population_caller.rs      Scaffold
      trio_caller.rs            Scaffold
      cancer_caller.rs          Scaffold
      polyclone_caller.rs       Scaffold
      cell_caller.rs            Scaffold
    models/
      haplotype_likelihood.rs   HaplotypeLikelihoodArray, BaseQualityLikelihoodModel,
                                FlatHaplotypeLikelihoodModel, phred_to_prob
      genotype/
        individual_model.rs     IndividualModel (log-sum-exp Bayesian inference)
        genotype_prior_model.rs UniformGenotypePriorModel,
                                HardyWeinbergGenotypePriorModel
    tools/
      vargen/
        cigar_scanner.rs        CigarScanner — generates candidate variants from pileup
        variant_generator.rs    VariantGenerator trait
      haplotype_filter.rs       Haplotype quality filter scaffold
      read_assigner.rs          Read→haplotype assignment scaffold
      ...
  io/
    reference/reference_genome.rs  ReferenceGenome + ReferenceReader trait
    variant/vcf_record.rs          VcfRecord, VcfRecordBuilder
    read/bam_reader.rs             BAM reader stub (htslib integration pending)
    region/region_parser.rs        Genomic region string parser
  readpipe/                    Read pipeline (filters, transformers, downsampling)
  containers/                  Mappable data structures
  Cargo.toml                   Package manifest
scripts/
  build_rust.sh                Smart build script (restores only missing files from git)
```

## IndividualCaller Pipeline (implemented)
1. **CigarScanner** scans all reads against reference → candidate SNVs and indels
2. **build_alt_haplotype()** splices each variant into the reference sequence
3. **HaplotypeLikelihoodArray** populated using `BaseQualityLikelihoodModel`
   (per-base Phred quality, genomic overlap, log(1−ε) match / log(ε/3) mismatch)
4. **make_genotypes()** enumerates all diploid combinations-with-repetition
5. **HardyWeinbergGenotypePriorModel** computes P(G) = C(n,k)·θ^k·(1−θ)^(n−k)
6. **IndividualModel.evaluate()** returns MAP genotype and posterior probabilities
7. Variant posterior = 1 − P(hom_ref); Phred-filtered at `min_variant_quality`
8. **VcfRecord** emitted per alt allele with GT / GQ / DP / AD FORMAT fields;
   anchor bases added for insertions/deletions per VCF specification

## CallerOptions defaults
| Field                  | Default  | Meaning                                |
|------------------------|----------|----------------------------------------|
| min_variant_quality    | 2.0      | Phred threshold for emitting a call    |
| max_haplotypes         | 128      | Cap on haplotype pool size             |
| min_read_depth         | 1        | Minimum reads to attempt a call        |
| ploidy                 | 2        | Sample ploidy                          |
| min_base_quality       | 20       | Phred threshold for candidate scanning |
| snp_heterozygosity     | 0.001    | θ for H-W SNP prior                    |
| indel_heterozygosity   | 0.0001   | θ for H-W indel prior (stored)         |

## Next Steps
- **ReadManager**: integrate `rust-htslib` for real BAM parsing
  (`PKG_CONFIG_PATH=/nix/store/29ardlwynaqws6ay3lvwdds5wjh3h4qr-htslib-1.15/lib/pkgconfig`)
- **Population/Trio/Cancer callers**: implement calling logic analogous to IndividualCaller
- **HMM likelihood model**: full pair-HMM alignment for indel accuracy
- **VCF writer**: write multi-sample VCF header + sorted records to output file

## Session Notes
- Files are in git but NOT on disk after a session restart.
  `scripts/build_rust.sh` restores only missing files (safe for uncommitted edits).
- Always set `export PATH="$HOME/.nix-profile/bin:$PATH"` before any cargo commands.
- VcfRecord uses `NucleotideSequence = String`; core types use `Vec<u8>`.
  Convert with `String::from_utf8_lossy(&bytes).to_string()` at the VCF layer.
