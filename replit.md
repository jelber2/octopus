# Octopus Variant Caller — Rust Port

## Project Overview

Full conversion of the **Octopus** bioinformatics variant caller from C++ to Rust.
- Original C++: `src/` (636 files, ~154,000 lines) — untouched
- Rust port: `rust_src/` — Cargo.toml at workspace root
- Build command: `cargo build 2>&1` (workflow: "Build Octopus (Rust)")

## Build Status

**Clean build** — zero errors, only expected "unused code" warnings.
Last verified: 2026-03-15

## CLI Argument Parsing (`rust_src/main.rs`)

`main.rs` implements the full Octopus CLI via `clap` (derive API), matching all
argument groups from the original `src/config/option_parser.cpp`.

| Group | Key flags |
|---|---|
| General | `-R/--reference` (required), `-I/--reads` or `-i/--reads-file` (one required), `-T/--regions`, `-t/--regions-file`, `-K/--skip-regions`, `-o/--output`, `-S/--samples`, `--pedigree`, `--threads` (`--tentacles`), `--fast`, `--very-fast`, `--debug`, `--trace` |
| Read preprocessing | `--min-mapping-quality` (5), `--good-base-quality` (20), `--min-good-bases` (20), `--max-read-length` (10000), `--downsample-above` (1000), `--downsample-target` (500) |
| Variant discovery | `--variant-discovery-mode` (illumina), `--min-pileup-base-quality` (20), `--max-variant-size` (2000), `--kmer-sizes` (10 15 20), `--source-candidates` |
| Variant calling | `-C/--caller` (population), `-P/--organism-ploidy` (2), `-z/--snp-heterozygosity` (0.001), `-y/--indel-heterozygosity` (0.0001), `--min-variant-posterior` (0.1), `--sequence-error-model`, `--min-phase-score` (5.0) |
| Cancer | `-N/--normal-samples`, `--somatic-snv-prior`, `--somatic-indel-prior`, `--min-credible-somatic-frequency`, `--somatics-only` |
| Trio | `-M/--maternal-sample`, `-F/--paternal-sample`, `--denovo-snv-prior`, `--denovo-indel-prior`, `--denovos-only` |
| Polyclone | `--max-clones` (5), `--min-clone-frequency`, `--clone-prior` |
| Cell | `--max-copy-loss`, `--max-copy-gain`, `--dropout-concentration` |

`--caller` is validated at parse time against the 6 model names.
`validate_args()` checks all input files exist and flag conflicts (`--fast` + `--very-fast`,
downsampler target ≥ above, etc.) and emits caller-specific errors (cancer requires
`--normal-samples`; trio requires both `--maternal-sample` and `--paternal-sample`).

## Architecture

### Module Structure
```
rust_src/
  main.rs               # Entry point, declares all modules
  basics/               # Fundamental types (GenomicRegion, AlignedRead, CigarString, etc.)
  concepts/             # Mappable, Comparable traits
  utils/                # Math, string utils, thread pool, memoize
  exceptions/           # OctopusError, Result type
  config/               # Configuration structs
  logging/              # Log levels and formatting
  io/
    reference/          # ReferenceGenome trait + FASTA reader
    variant/            # VcfRecord, VcfHeader, VcfReader, VcfWriter, VcfSpec
    read/               # BamReader stub, ReadManager
    region_parser.rs    # Region string parsing
  core/
    types/              # Allele, Variant, Haplotype, Genotype, CancerGenotype, Phylogeny
    models/
      haplotype_likelihood/  # PairHMM, FlatHaplotypeLikelihoodModel, HaplotypeLikelihoodArray
      genotype/              # Individual, Population, Trio, Cancer, Subclone, HardyWeinberg models
      mutation/              # Somatic and de-novo mutation models
    callers/            # 6 callers: Individual, Population, Trio, Cancer, Polyclone, Cell
    csr/                # Call set refinement: facets, filters, measures
    tools/
      hapgen/           # HaplotypeGenerator, HaplotypeTree, GenomeWalker
      vargen/           # CigarScanner, LocalReassembler, VariantGenerator trait
      phaser/           # Phaser, PhaseSet
      vcf_header_factory.rs
      vcf_record_factory.rs
      haplotype_filter.rs
      read_assigner.rs
      genome_walker.rs
      bam_realigner.rs
      bad_region_detector.rs
  containers/           # MappableFlatSet, MappableFlatMultiSet, MappableMap, MappableBlock, ProbabilityMatrix
  readpipe/
    filtering/          # ReadFilter trait + common filters (MappingQuality, ReadLength, etc.)
    transformers/       # ReadTransformer trait + common transformers
    downsampling/       # Downsampler trait + CoverageDownsampler
    read_pipe.rs        # ReadPipe orchestrator
```

### Key Design Decisions
- `boost::optional` → `Option<T>`
- C++ exceptions → `Result<T, OctopusError>`
- C++ templates → Rust generics + traits
- `ReferenceGenome` is a trait for polymorphic reference backends (FASTA, etc.)
- `NucleotideSequence = Vec<u8>` in core types; `String` in VCF record types
- `CigarString = Vec<CigarOperation>` with `CigarFlag` enum for CIGAR operations
- `MappableFlatSet::overlap_range` uses explicit `'a: 'b` lifetime bounds for iterator capture
- `ProbabilityMatrix<R, C>` uses `PhantomData<(R, C)>` to satisfy unused type parameter constraints

## Dependencies (Cargo.toml)
- `log` + `env_logger` — logging
- Standard library only otherwise (no heavy bio crates)
