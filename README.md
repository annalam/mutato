# Mutato

Mutato is a software package for detecting somatic mutations and germline variants in DNA sequencing data. It takes as input a number of samples in BAM format, and outputs a table listing all genomic positions where variant alleles were found, and the number of variant allele reads in every sample. The availability of variant allele read count information from all samples allows easy implementation of mutation calling pipelines where known negative control samples are used to generate an allele-specific background error rate for each mutation.

Features
--------
- Simultaneous variant calling across multiple BAM files
- High performance (analyzes xxx,xxx gigabytes of BAM format alignments / hour / thread)
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
mutato call --alt-reads=5 --alt-frac=0.05 --min-mapq=0 hg38.fa *.bam > variants.vcf
```

The first positional parameter `hg38.fa` is a reference genome FASTA file, followed by a list of BAM files to analyze. The BAM files must be position-sorted, and must have been aligned against the provided reference genome.  All samples of interest (including tumor samples, matched germline samples, and negative control samples) should be included in the same run. Mutato will detect all candidate variants that carry a non-reference allele with read level evidence surpassing the user-specified thresholds `--alt-reads` and `--alt-frac`. The generated output file is a tab-delimited matrix describing the read-level evidence for each candidate variant in each input BAM file.

The analysis can be parallelized across chromosomes by utilizing the `--region` flag:
```
mutato call --region=chr1 --alt-reads=5 --alt-frac=0.05 --min-mapq=0 hg38.fa *.bam > variants.vcf
```


Future work
-----------
- Add functionality for post-processing and annotating somatic mutations and germline variants into the repository (currently implemented as Julia language scripts in our own pipelines)

