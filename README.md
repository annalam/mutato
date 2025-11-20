# Mutato

Mutato is a software package for detecting somatic mutations and germline variants in DNA sequencing data. It takes as input a number of samples in BAM or CRAM format, and outputs a TSV-format table listing all genomic positions where non-reference alleles were found, as well as numerical data and statistics about the evidence supporting those variants. The availability of variant allele read count information from all samples allows easy implementation of mutation calling pipelines where known negative control samples are used to generate an allele-specific background error rate for each mutation.

Features
--------
- Multi-threaded parallel variant calling across multiple BAM/CRAM samples
- Allows empirical background error estimation based on read-level evidence from a set of background control samples
- Implemented entirely in the Rust language (no unsafe code)
- Self-contained binary, no external dependencies


Installation
------------

Install Rust (version 1.44 or later). Then run the following command:
```
cargo install --git https://github.com/annalam/mutato
```

Usage
-----

First run `mutato call` on a set of input BAM files to detect candidate variants and to generate a matrix of read-level evidence:
```
mutato call --alt-reads=5 --alt-frac=0.05 hg38.fa *.bam > variants.vcf
```

The first positional parameter `hg38.fa` is a reference genome FASTA file, followed by a list of BAM/CRAM files to analyze. The BAM/CRAM files must be position-sorted.  All samples of interest (including tumor samples, matched germline samples, and negative control samples) are included in the same run. Mutato will detect all candidate variants that carry a non-reference allele with read level evidence surpassing the user-specified thresholds `--alt-reads` and `--alt-frac`. The generated output file is a tab-delimited matrix describing the read-level evidence for each candidate variant in each input BAM/CRAM file.
