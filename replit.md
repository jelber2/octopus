# Octopus Variant Caller — Rust Port

## Project Overview

Full conversion of the **Octopus** bioinformatics variant caller from C++ to Rust.
- Original C++: `src/` (636 files, ~154,000 lines) — untouched
- Rust port: `rust_src/` — Cargo.toml at workspace root
- Build command: `cargo build 2>&1` (workflow: "Build Octopus (Rust)")

## Build Status

**Clean build** — zero errors, only expected "unused code" warnings (~461).
Last verified: 2026-03-15

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
