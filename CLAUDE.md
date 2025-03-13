# CLAUDE.md - Coding Assistant Guide

## Core Concepts
- Genomic data analysis pipeline built with Snakemake workflow manager
- Processes sequencing data through extraction, alignment, and feature overlap analysis
- DAG (Directed Acyclic Graph) workflow with multiple dependent stages

## Build Commands
- Run pipeline: `uv run snakemake -c1`
- Clean outputs: `uv run snakemake clean`
- Run specific rule: `uv run snakemake -c1 <rule_name>`
- Run Python script: `uv run python <script_name.py>`
- Run Jupyter notebook: `uv run jupyter notebook report.ipynb`

## Dependencies
- Package management via `uv` (see `pyproject.toml` and `uv.lock`)
- Command-line tools required:
  - nanoq (for quality control)
  - cutadapt (for adapter trimming)
  - seqtk (for sequence format conversion)
  - minimap2 (for sequence alignment)
  - samtools (for BAM file processing)
  - bedtools (for genomic feature overlap)

## Code Style Guidelines
- **Imports**: Group standard, third-party, and local imports
- **Functions**: Well-documented with docstrings
- **Error handling**: Use try/finally for cleanup operations
- **Output paths**: Consistent directory structure in `output/`
- **Configuration**: Parameters stored in `config.yaml`