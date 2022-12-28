use std::str;
use rayon::prelude::*;
use crate::common::{parse_args, Genome, FileReader};
use crate::call::{Mutation, Settings, analyze_bam_secondpass};

const USAGE: &str = "
Usage:
  mutato call targeted [options] <genome.fa> <variants.tsv> <bam_files>...

Options:
  --min-baseq=N     Ignore bases with a base quality < N [default: 0]
  --min-end-dist=N  Ignore N bases at both read ends [default: 0]
  --no-strand-bias  Do not report strand bias statistics
  --threads=N       Maximum number of threads to use [default: 4]

Calculates the number of supporting reads for a list of genomic variants,
in a set of BAM files provided by the user. The following additional statistics
are reported for each variant:
- Average mapping quality (MAPQ) of reads supporting the variant allele
- Average distance of the variant allele from the nearest read end
- Percentage of variant allele reads aligned to the + strand of the genome
- Percentage of non-variant allele reads aligned to the + strand of the genome
";

pub fn main() {
	let args = parse_args(USAGE);
	let genome_path = args.get_str("<genome.fa>");
	let variants_path = args.get_str("<variants.tsv>");
	let bam_paths = args.get_vec("<bam_files>");
	let max_threads: usize = args.get_str("--threads").parse().unwrap_or_else(
		|_| error!("--threads must be a positive integer."));

	let mut settings = Settings::default();
	settings.alt_reads = 0;
	settings.alt_frac = 0.0;
	settings.min_baseq = args.get_str("--min-baseq").parse().unwrap_or_else(
		|_| error!("--min-baseq must be a non-negative integer."));
	settings.min_end_distance = args.get_str("--min-end-dist").parse()
		.unwrap_or_else(|_| error!("--min-end-dist must be a non-negative integer."));
	settings.no_strand_bias = args.get_bool("--no-strand-bias");

	eprintln!("Reading reference genome into memory...");
	let genome = Genome::from_fasta(&genome_path);

	// Initialize the thread pool for analyzing multiple BAM files in parallel
	rayon::ThreadPoolBuilder::new().num_threads(max_threads).build_global()
		.unwrap();

	eprintln!("Reading the variant list into memory...");
	let mut variants_tsv = FileReader::new(&variants_path);
	let mut variants: Vec<Mutation> = Vec::new();

	let mut line = String::new();
	while variants_tsv.read_line(&mut line) {
		if line.starts_with('#') { continue; }
		let cols: Vec<&str> = line.trim_end().split('\t').collect();
		if cols.len() < 4 { continue; }
		let chr = genome.chr_idx(&cols[0]);
		let position: u32 = cols[1].parse().unwrap_or_else(|_|
			error!("Variant position '{}' is invalid:\n{}", &cols[1], &line));

		// We need to convert the alternate allele into the internal
		// representation used by Mutato.
		variants.push(Mutation { chr: chr as u16, position, alt_allele:
			if cols[3].len() == 1 && cols[2].len() == 1 {
				cols[3].into()
			} else if cols[3].len() > 1 && cols[2].len() == 1 {
				format!("+{}", &cols[3][1..]).into()
			} else if cols[2].len() > 1 && cols[3].len() == 1 {
				let mut del: String = "-".into();
				for k in 1..cols[2].len() { del.push('N'); }
				del.into()
			} else {
				eprintln!("ERROR: '{}' > '{}'", &cols[2], &cols[3]);
				unreachable!()
			}
		});
	}

	eprintln!("Calculating evidence matrix for {} variants...", variants.len());
	let evidence: Vec<_> = bam_paths.par_iter().map(|&bam_path| {
		analyze_bam_secondpass(&bam_path, &genome, &variants, &settings)
	}).collect();

	eprintln!("Analysis complete.");

	// This code is duplicated from call.rs.
	print!("CHROM\tPOSITION\tREF\tALT\tNOTES");
	for bam_path in &bam_paths {
		print!("\t{}", bam_path.trim_end_matches(".bam"));
	}
	println!();
	for k in 0..variants.len() {
		let chr = variants[k].chr;
		let pos = variants[k].position as usize;
		print!("{}\t{}\t", &genome.name(chr as usize), pos);

		let ref_base = genome.sequence_by_chr_idx(chr as usize)[pos - 1] as char;

		if variants[k].alt_allele.len() == 1 {
			print!("{}\t{}\t", ref_base, &variants[k].alt_allele);
		} else if variants[k].alt_allele.starts_with('+') {
			print!("{}\t{}{}\t", ref_base, ref_base,
				&variants[k].alt_allele[1..]);
		} else if variants[k].alt_allele.starts_with('-') {
			let del_len = variants[k].alt_allele.len() - 1;
			print!("{}\t{}\t", str::from_utf8(&genome.sequence_by_chr_idx(chr as usize)[pos-1..pos+del_len]).unwrap(), ref_base);
		} else {
			unreachable!();
		}

		for s in 0..bam_paths.len() {
			let e = &evidence[s][k];
			print!("\t{}:{}", e.alt_reads, e.total_reads);
			if e.alt_reads == 0 {
				print!(":0::0");
			} else {
				if settings.no_strand_bias {
					print!(":{}::{}", e.mapq, e.sidedness);
				} else {
					print!(":{}:{},{}:{}", e.mapq, e.strand, e.other_strand, e.sidedness);
				}
			}
		}
		println!();
	}
}

