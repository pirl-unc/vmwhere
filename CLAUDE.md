# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

VMwhere is a Python bioinformatics tool for discovering and genotyping tandem repeat (microsatellite) regions from long-read sequencing data. It has three commands: `find` (discover microsatellites in a reference FASTA), `genotype` (call alleles from BAM reads), and `visualize` (R-based allele frequency plots).

## Build and Development

```bash
# Install in editable mode with dev dependencies
pip install -e ".[dev]"

# Run all tests
pytest tests/ -v --tb=short

# Run a single test file
pytest tests/test_find.py -v

# Run a single test class or method
pytest tests/test_find.py::TestRunFind -v
pytest tests/test_find.py::TestRunFind::test_basic_find -v
```

Python >=3.12.4 required. R >=4.0 required only for the `visualize` subcommand.

## Architecture

### Module Layout (`src/vmwhere/`)

- **`cli.py`** — Argparse-based CLI dispatcher. Entry point registered as `vmwhere = vmwhere.cli:main`. Three subcommands: `find`, `genotype`, `visualize`.
- **`find.py`** — Microsatellite discovery pipeline. Scans a FASTA for exact motif matches (forward + reverse complement), counts consecutive repeats, filters by minimum repeat count, merges adjacent regions within a gap, and outputs a BED file.
- **`genotyper.py`** — Read extraction and genotyping. Fetches overlapping reads from BAM files, decomposes sequences using dynamic programming (edit-distance minimization against the motif), clusters reads by Levenshtein distance, and calls alleles based on read support thresholds. Uses `multiprocessing` to parallelize across chromosomes/motifs.
- **`visualize_region.R`** — R script invoked via subprocess by the `visualize` subcommand.

### Key Algorithm: Sequence Decomposition (`genotyper.py`)

The core genotyping algorithm decomposes each read into motif vs. non-motif segments using DP-based edit distance minimization. Key functions in the pipeline:
1. `decompose_seq_with_motif_anchors()` — splits sequence at exact motif matches
2. `decompose_non_motif_region()` / `further_parse_intervening()` — minimizes edit distance for intervening regions
3. `cluster_consolidated_reads_by_edit_distance()` — groups similar decomposed reads
4. `identify_alleles()` — calls major/minor alleles from clustered read support

### Data Flow

`find` produces a BED file → `genotype` reads that BED + BAM + FASTA → outputs a TSV with VCF-like columns (CHROM, POS, ID, REF, ALT, GT, CN, etc.) → `visualize` reads the TSV to produce PDF plots.

### Chromosome Name Mapping

`chr_mapping_simple.txt` (bundled as package data) maps RefSeq accession IDs to UCSC chromosome names (e.g., `NC_060925.1` → `chr1`). Used in `find.py` when processing T2T or RefSeq-style FASTA files.

## Testing

Tests are in `tests/` using pytest. `conftest.py` adds `src/` to `sys.path`. Test data files (BAM, BED, FASTA) live alongside tests. The main test file is `test_find.py` (~550 lines) covering unit and integration tests for the find pipeline. `test_decompose_reads.py` and `test_cluster_reads.py` cover genotyper functions.

## Versioning

Auto-versioned via `setuptools-scm` from git tags. Version is exported from `src/vmwhere/__init__.py` which reads `_version.py` (auto-generated, gitignored).

## CI

GitHub Actions runs pytest on Python 3.12 and 3.13 on push/PR to main. Releases are triggered by `v*` tags and publish to PyPI via trusted publisher (OIDC).
