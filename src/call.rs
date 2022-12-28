use std::cmp::min;
use std::str;
use std::collections::VecDeque;
use rust_htslib::bam::{self, Read, record::Cigar};
use itertools::Itertools;
use rayon::prelude::*;
use smartstring::alias::String;    // Small string optimization
use crate::common::{parse_args, Genome};

const USAGE: &str = "
Usage:
  mutato call [options] <genome.fa> <bam_files>...

Options:
  --alt-reads=N     Minimum alt allele reads to report [default: 5]
  --alt-frac=N      Minimum alt allele fraction to report [default: 0.1]
  --min-baseq=N     Ignore bases with a base quality < N [default: 0]
  --min-end-dist=N  Ignore N bases at both read ends [default: 0]
  --ffpe            Suppress formaldehyde induced artifacts (e.g. FFPE sample)
  --oxidative       Suppress oxidative damage (i.e. oxo-guanine)
  --no-strand-bias  Do not report strand bias statistics
  --threads=N       Maximum number of threads to use [default: 4]
  --single-pass     Generate candidate variant list, but skip the second pass

Analyzes a set of BAM samples provided by the user to identify potential
genomic variants. The analysis is carried out in two passes. In the first pass,
the algorithm searches for genomic positions where at least one of the BAM
samples carries a non-reference allele with the number of supporting reads
exceeding user-specified thresholds. In the second pass, the algorithm
calculates the precise number of (non-duplicate) reads supporting each
candidate variant in each BAM file, and outputs this information as a
tab-separated file reminiscent of the VCF format.

The software also reports the following additional statistics for each
candidate variant:
- Average mapping quality (MAPQ) of reads supporting the variant allele
- Average distance of the variant allele from the nearest read end
- Percentage of variant allele reads aligned to the + strand of the genome
- Percentage of non-variant allele reads aligned to the + strand of the genome

By including negative control samples in the analysis, an empirical background
error rate can be established for each variant, and used for filtering.
";


#[derive(Copy, Clone)]
struct Substitution {
	reads: u32,
	sidedness: u32,     // Sum of distances from nearest read end
	mapq: u32,          // Sum of MAPQ values of reads supporting this allele
	read_strand: u32,   // Supporting reads in + strand
	phys_strand: u32    // Supporting reads from ssDNA fragments in + strand
}

impl Substitution {
	pub fn new() -> Substitution {
		Substitution { reads: 0, sidedness: 0, mapq: 0, read_strand: 0, phys_strand: 0 }
	}
}

#[derive(Clone)]
struct Indel {
	sequence: String,   // Prefixed with '+' for insertion, '-' for deletion
	reads: u32,
	sidedness: u32,     // Average position from nearest read end
	mapq: u32,
	read_strand: u32,
	phys_strand: u32
}

#[derive(Clone)]
struct Pileup {
	total: u32,
	acgt: [Substitution; 4],
	indels: Vec<Indel>,
	read_strand: u32,       // Total number of overlapping reads in + strand
	phys_strand: u32
}

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct Mutation {
	pub chr: u16,
	pub position: u32,
	pub alt_allele: String
}

#[derive(Clone)]
pub struct Evidence {
	pub alt_reads: u32,
	pub total_reads: u32,
	pub mapq: u8,
	pub sidedness: u8,
	pub strand: u8,         // Percent of alt allele reads in + strand (0..100)
	pub other_strand: u8    // Percent of other reads in + strand (0..100)
}

#[derive(Default)]
pub struct Settings {
	pub alt_reads: u32,
	pub alt_frac: f32,
	pub min_baseq: u8,
	pub min_end_distance: u8,
	pub no_strand_bias: bool,
	pub suppress_ffpe_artifacts: bool,
	pub suppress_oxidative_damage: bool
}

pub fn main() {
	let args = parse_args(USAGE);
	let genome_path = args.get_str("<genome.fa>");
	let bam_paths = args.get_vec("<bam_files>");
	let max_threads: usize = args.get_str("--threads").parse().unwrap_or_else(
		|_| error!("--threads must be a positive integer."));
	let single_pass = args.get_bool("--single-pass");

	let mut settings = Settings::default();
	settings.alt_reads = args.get_str("--alt-reads").parse().unwrap_or_else(
		|_| error!("--alt-reads must be a positive integer"));
	settings.alt_frac = args.get_str("--alt-frac").parse().unwrap_or_else(
		|_| error!("--alt-frac must be a fraction between 0 and 1."));
	settings.min_baseq = args.get_str("--min-baseq").parse().unwrap_or_else(
		|_| error!("--min-baseq must be a non-negative integer."));
	settings.min_end_distance = args.get_str("--min-end-dist").parse()
		.unwrap_or_else(|_| error!("--min-end-dist must be a non-negative integer."));
	
	settings.no_strand_bias = args.get_bool("--no-strand-bias");
	settings.suppress_ffpe_artifacts = args.get_bool("--ffpe");
	settings.suppress_oxidative_damage = args.get_bool("--oxidative");

	eprintln!("Reading reference genome into memory...");
	let genome = Genome::from_fasta(&genome_path);

	// Initialize the thread pool for analyzing multiple BAM files in parallel
	rayon::ThreadPoolBuilder::new().num_threads(max_threads).build_global()
		.unwrap();

	eprintln!("Analyzing BAM files for candidate variants:");
	let mut sample_variants: Vec<_> = bam_paths.par_iter().map(|&bam_path| {
		eprintln!("- {}", &bam_path);
		let mut variants = analyze_bam_firstpass(&bam_path, &genome, &settings);
		variants.sort();  // Sort in preparation for k-way merge
		variants
	}).collect();

	eprintln!("Merging per-sample lists of candidate variants...");
	let variants: Vec<Mutation> = itertools::kmerge(sample_variants).dedup().collect();

	eprintln!("Found {} candidate variants.", variants.len());

	eprintln!("Calculating evidence matrix...");

	let evidence: Vec<_> = bam_paths.par_iter().map(|&bam_path| {
		if single_pass { return Vec::new(); }  // Skip the second pass
		analyze_bam_secondpass(&bam_path, &genome, &variants, &settings)
	}).collect();
	
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

		if single_pass { println!(); continue; }   // Skip the second pass

		for s in 0..evidence.len() {
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


fn analyze_bam_firstpass(bam_path: &str, genome: &Genome, settings: &Settings) -> Vec<Mutation> {

	// Open the BAM file for reading
	let mut bam = bam::Reader::from_path(bam_path).unwrap_or_else(
		|_| error!("Could not open BAM file '{}'.", &bam_path));
	let header = bam.header().clone();
	let chr_names: Vec<&str> = header.target_names().iter().map(|x| str::from_utf8(x).unwrap()).collect();

	let mut pileups: VecDeque<Pileup> = VecDeque::with_capacity(1000);
	let mut curr_chr: u32 = u32::MAX;
	let mut curr_pos: u32 = 0;
	let mut chr_seq = genome.sequence_by_chr_idx(0);
	let mut mutations: Vec<Mutation> = Vec::new();

	let mut read = bam::Record::new();
	while read_bam_record(&mut bam, &mut read) {
		if read.is_unmapped() || read.is_duplicate() { continue; }
		if read.is_secondary() || read.is_supplementary() { continue; }
		if read.is_quality_check_failed() { continue; }
		if read.tid() < 0 { error!("Invalid TID < 0."); }

		let chr = read.tid() as u32;
		let pos = read.pos() as u32 + 1;   // Convert to one-based coordinate

		// Check how many pileups we have finalized (i.e. will no longer change)
		let to_report = if chr != curr_chr {
			pileups.len()
		} else if pos > curr_pos {
			min((pos - curr_pos) as usize, pileups.len())
		} else if pos < curr_pos {
			error!("BAM file must be sorted by position.");
		} else {
			0
		};

		// Go through the finalized pileups, and check if any meet our
		// thresholds.
		for pileup in pileups.drain(0..to_report) {
			let ref_base = chr_seq[curr_pos as usize - 1];
			let chr_idx = genome.chr_idx(&chr_names[curr_chr as usize]) as u16;

			const ACGT: [u8; 4] = [b'A', b'C', b'G', b'T'];
			for base in 0..4 {
				if ref_base == ACGT[base] { continue; }
				let reads = pileup.acgt[base].reads;
				if reads < settings.alt_reads || (reads as f32) / (pileup.total as f32) < settings.alt_frac {
					continue;
				}
				mutations.push(Mutation { chr: chr_idx, position: curr_pos,
					alt_allele: str::from_utf8(&[ACGT[base]]).unwrap().into() });
			}

			for indel in pileup.indels {
				if indel.reads < settings.alt_reads || (indel.reads as f32) /
					(pileup.total as f32) < settings.alt_frac {
					continue;
				}
				let sign = indel.sequence.as_bytes()[0];
				mutations.push(Mutation { chr: chr_idx, position: curr_pos,
					alt_allele: indel.sequence.into() });
			}

			curr_pos += 1;
		}

		// It is critical that these are updated only *after* reporting
		// the pileups.
		if chr != curr_chr {
			chr_seq = genome.sequence(chr_names[chr as usize]);
		}
		curr_chr = chr;
		curr_pos = pos;

		// At this point pileups[0] represents chromosome position "curr_pos".
		// We add the read to the pileups vector.
		add_read_to_pileups(&mut pileups, &read, &settings);
	}

	mutations
}



pub fn analyze_bam_secondpass(bam_path: &str, genome: &Genome, mutations: &Vec<Mutation>, settings: &Settings) -> Vec<Evidence> {

	// Open the BAM file for reading
	let mut bam = bam::Reader::from_path(bam_path).unwrap_or_else(
		|_| error!("Could not open BAM file '{}'.", &bam_path));
	let header = bam.header().clone();

	// For each chromosome listed in the BAM header, calculate its index
	// within the Genome data structure.
	let chr_indices: Vec<u16> = header.target_names().iter().map(|x| genome.chr_idx(&str::from_utf8(x).unwrap()) as u16).collect();

	let mut pileups: VecDeque<Pileup> = VecDeque::with_capacity(1000);
	let mut curr_chr: u32 = u32::MAX;
	let mut curr_pos: u32 = 0;
	let mut chr_seq = genome.sequence_by_chr_idx(0);
	let mut chr_muts: VecDeque<usize> = VecDeque::new();
	let mut evidence = vec![Evidence { alt_reads: 0, total_reads: 0, mapq: 0, sidedness: 0, strand: 0, other_strand: 0 }; mutations.len()];

	let mut read = bam::Record::new();
	while read_bam_record(&mut bam, &mut read) {
		if read.is_unmapped() || read.is_duplicate() { continue; }
		if read.is_secondary() || read.is_supplementary() { continue; }
		if read.is_quality_check_failed() { continue; }
		if read.tid() < 0 { error!("Invalid TID < 0."); }

		let chr = read.tid() as u32;
		let pos = read.pos() as u32 + 1;   // Convert to 1-based coordinate

		// Check how many pileups we have finalized (i.e. will no longer change)
		let to_report = if chr != curr_chr {
			pileups.len()
		} else if pos > curr_pos {
			min((pos - curr_pos) as usize, pileups.len())
		} else if pos < curr_pos {
			error!("BAM file must be sorted by position.");
		} else {
			0
		};

		// Go through the finalized pileups, and check if the user requested
		// any of them to be reported out.
		while !chr_muts.is_empty() &&
			mutations[chr_muts[0]].position < curr_pos {
			chr_muts.pop_front();
		}
		for pileup in pileups.drain(0..to_report) {
			while !chr_muts.is_empty() &&
				mutations[chr_muts[0]].position == curr_pos {

				// This genomic position contained a candidate variant in the
				// first-pass analysis, so we store the number of supporting
				// reads in the matrix.
				let mutation = &mutations[chr_muts[0]];
				evidence[chr_muts[0]] = calculate_mutation_evidence(
					&mutation, &pileup, &genome,
					settings.suppress_ffpe_artifacts,
					settings.suppress_oxidative_damage);
				chr_muts.pop_front();
			}
			curr_pos += 1;
		}

		// If we have moved to a new chromosome, we generate a new vector
		// containing the indices of that chromosome's mutation candidates,
		// sorted in ascending order of chromosomal position.
		if chr != curr_chr {
			let mut inside: Vec<_> = (0..mutations.len()).filter(|&m| mutations[m].chr == chr_indices[chr as usize]).collect();
			inside.sort_by_key(|&m| mutations[m].position);
			chr_muts = inside.into();   // Convert from Vec to VecDeque

			chr_seq = genome.sequence_by_chr_idx(chr_indices[chr as usize] as usize);
		}

		// It is critical that these are updated only *after* reporting
		// the pileups.
		curr_chr = chr;
		curr_pos = pos;

		// At this point pileups[0] represents chromosome position "curr_pos".
		// We add the read to the pileups vector.
		add_read_to_pileups(&mut pileups, &read, &settings);
	}

	evidence
}



fn add_read_to_pileups(pileups: &mut VecDeque<Pileup>, read: &bam::Record, settings: &Settings) {
	// Calculate how long the pileup vector needs to be to accommodate
	// the information from this read.
	let mut read_span = 0;
	for s in read.cigar().iter() {
		match *s {
			Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) |
			Cigar::Del(len) | Cigar::RefSkip(len) => {
				read_span += len as usize;
			},
			Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) => {},
			_ => error!("Unsupported CIGAR element found.")
		}
	}

	// Extend the pileup vector if necessary.
	while pileups.len() < read_span {
		pileups.push_back(Pileup {
			acgt: [Substitution { reads: 0, mapq: 0, sidedness: 0, read_strand: 0, phys_strand: 0 }; 4],
			total: 0, indels: Vec::new(), read_strand: 0, phys_strand: 0 });
	}

	// These variables keep track of our current position
	// within the pileup track and the read sequence, as we
	// read through the CIGAR string.
	let mut pileup_idx = 0;
	let mut seq_idx = 0;
	let mut prior_n = false;

	let seq = read.seq();
	let qual = read.qual();

	// In typical Illumina paired end sequencing + library preparation,
	// the mate #1 strand represents the original physical genome strand of
	// the single-stranded DNA fragment that was sequenced. Here we figure
	// out that physical strand.
	let native_minus = read.is_reverse() == read.is_first_in_template();

	for s in read.cigar().iter() {
		match *s {
			Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
				for _ in 0..len {
					let acgt_idx = encoded_base_to_acgt_index(seq.encoded_base(seq_idx));
					if acgt_idx >= 4 {
						prior_n = true;    // Ambiguous nucleotide
					} else if qual[seq_idx] < settings.min_baseq {
						prior_n = true;    // Low quality base, ignore
					} else {
						let mut pileup = &mut pileups[pileup_idx];
						let mut sub = &mut pileup.acgt[acgt_idx];

						prior_n = false;

						let end_distance = min(seq_idx, seq.len() - 1 - seq_idx) as u32;
						if end_distance >= settings.min_end_distance as u32 {
							pileup.total += 1;
							pileup.read_strand += if read.is_reverse() { 0 } else { 1 };
							pileup.phys_strand += if native_minus { 0 } else { 1 };

							sub.reads += 1;
							sub.read_strand += if read.is_reverse() { 0 } else { 1 };
							sub.phys_strand += if native_minus { 0 } else { 1 };

							// Running totals, division later
							sub.sidedness += end_distance;
							sub.mapq += read.mapq() as u32;
						}
					}
					pileup_idx += 1;
					seq_idx += 1;
				}
			},
			Cigar::Ins(len) => {
				// TODO: Add support for insertions at the beginning of a read.
				// Currently ignored to prevent indexing pileups[-1].
				if pileup_idx == 0 { seq_idx += len as usize; continue; }

				let mut allele = String::new();
				allele.push('+');     // Mark insertions with a plus
				for _ in 0..len {
					allele.push(encoded_base_to_char(seq.encoded_base(seq_idx)));
					seq_idx += 1;
				}
				if allele.contains('N') { continue; }  // No ambiguous

				// Distance from nearest read end
				let end_distance = min(seq_idx - (len as usize),
					seq.len() - 1 - seq_idx) as u32;
				if end_distance < settings.min_end_distance as u32 {
					continue;
				}

				let pileup = &mut pileups[pileup_idx - 1];
				count_indel(pileup, &allele, end_distance, read.mapq(),
					!read.is_reverse(), !native_minus);

				// Indel read counts are stored in the pileup structure that
				// is located immediately to the left of where the indel begins.
				// If that base was N in this read, it was not counted towards
				// pileup.total and so we do so here.
				if prior_n {
					pileup.total += 1;
					pileup.read_strand += if read.is_reverse() { 0 } else { 1 };
					pileup.phys_strand += if native_minus { 0 } else { 1 };
				}
			},
			Cigar::Del(len) => {
				if pileup_idx == 0 { 
					error!("CIGAR strings starting with deletion are not supported.");
				}

				// Distance from nearest read end
				let mut pileup = &mut pileups[pileup_idx - 1];

				let mut allele = String::new();
				allele.push('-');
				for _ in 0..len {
					allele.push('N');
					pileup_idx += 1;
				}

				let end_distance = min(seq_idx, seq.len() - 1 - seq_idx) as u32;
				if end_distance < settings.min_end_distance as u32 {
					continue;
				}

				count_indel(&mut pileup, &allele, end_distance, read.mapq(),
					!read.is_reverse(), !native_minus);

				// Indel read counts are stored in the pileup structure that
				// is located immediately to the left of where the indel begins.
				// If that base was N in this read, it was not counted towards
				// pileup.total and so we do so here.
				if prior_n {
					pileup.total += 1;
					pileup.read_strand += if read.is_reverse() { 0 } else { 1 };
					pileup.phys_strand += if native_minus { 0 } else { 1 };
				}
			},
			Cigar::RefSkip(len) => { pileup_idx += len as usize; },
			Cigar::SoftClip(len) => { seq_idx += len as usize; },
			Cigar::HardClip(_) => {},
			_ => error!("Unsupported CIGAR element")
		}
	}
}




fn count_indel(pileup: &mut Pileup, seq: &str, end_distance: u32, mapq: u8, plus_read_strand: bool, plus_phys_strand: bool) {
	for i in 0..pileup.indels.len() {
		let mut indel = &mut pileup.indels[i];
		if &*indel.sequence == seq {
			indel.reads += 1;
			indel.sidedness += end_distance;
			indel.mapq += mapq as u32;
			indel.read_strand += if plus_read_strand { 1 } else { 0 };
			indel.phys_strand += if plus_phys_strand { 1 } else { 0 };
			return;
		}
	}
	pileup.indels.push(Indel {
		sequence: seq.into(), 
		reads: 1, 
		sidedness: end_distance as u32, 
		mapq: mapq as u32,
		read_strand: if plus_read_strand { 1 } else { 0 },
		phys_strand: if plus_phys_strand { 1 } else { 0 }
	});
}


fn calculate_mutation_evidence(mutation: &Mutation, pileup: &Pileup,
	genome: &Genome, fix_ffpe: bool, fix_oxog: bool) -> Evidence {
	
	if mutation.alt_allele.len() == 1 {
		// Base substitution
		let alt_base: u8 = mutation.alt_allele.as_bytes()[0];
		let base_idx = base_to_acgt_index(mutation.alt_allele.as_bytes()[0]);
		if base_idx == 10 {
			error!("Invalid alt allele '{}'.", mutation.alt_allele);
		}
		let sub = &pileup.acgt[base_idx];
		let mut evidence = Evidence {
			alt_reads: sub.reads, total_reads: pileup.total,
			mapq: (sub.mapq as f32 / sub.reads as f32).round() as u8,
			sidedness: (sub.sidedness as f32 / sub.reads as f32).round() as u8,
			strand: (sub.read_strand as f32 / sub.reads as f32 * 100.0).round() as u8,
			other_strand: ((pileup.read_strand - sub.read_strand) as f32 / (pileup.total - sub.reads) as f32 * 100.0).round() as u8
		};

		let ref_base = genome.sequence_by_chr_idx(mutation.chr as usize)[mutation.position as usize - 1];
		
		// Suppress formaldehyde-induced DNA damage (C>T) and
		// oxidative DNA damage (G>T) by only counting sequenced ssDNA
		// fragments that originate from the strand that cannot be affected
		// by these forms of single-stranded damage.
		if fix_ffpe && ref_base == b'C' && alt_base == b'T' {
			// Only count ssDNA fragments representing minus strand
			evidence.alt_reads = sub.reads - sub.phys_strand;
			evidence.total_reads = pileup.total - pileup.phys_strand;
		} else if fix_ffpe && ref_base == b'G' && alt_base == b'A' {
			// Only count ssDNA fragments representing plus strand
			evidence.alt_reads = sub.phys_strand;
			evidence.total_reads = pileup.phys_strand;
		} else if fix_oxog && ref_base == b'G' && alt_base == b'T' {
			// Only count ssDNA fragments representing minus strand
			evidence.alt_reads = sub.reads - sub.phys_strand;
			evidence.total_reads = pileup.total - pileup.phys_strand;
		} else if fix_oxog && ref_base == b'C' && alt_base == b'A' {
			// Only count ssDNA fragments representing plus strand
			evidence.alt_reads = sub.phys_strand;
			evidence.total_reads = pileup.phys_strand;
		}

		evidence

	} else if let Some(indel) = pileup.indels.iter()
		.find(|i| i.sequence == mutation.alt_allele) {
		Evidence {
			alt_reads: indel.reads, total_reads: pileup.total,
			mapq: (indel.mapq as f32 / indel.reads as f32).round() as u8,
			sidedness: (indel.sidedness as f32 / indel.reads as f32).round() as u8,
			strand: (indel.read_strand as f32 / indel.reads as f32 * 100.0).round() as u8,
			other_strand: ((pileup.read_strand - indel.read_strand) as f32 / (pileup.total - indel.reads) as f32 * 100.0).round() as u8
		}
	} else {
		Evidence { alt_reads: 0, total_reads: pileup.total,
			mapq: 0, sidedness: 0, strand: 0, other_strand: 0 }
	}
}






fn encoded_base_to_acgt_index(encoded: u8) -> usize {
	match encoded {
		1 => 0,   // A
		2 => 1,   // C
		4 => 2,   // G
		8 => 3,   // T
		_ => 4    // N
	}
}


fn encoded_base_to_char(encoded: u8) -> char {
	match encoded { 1 => 'A', 2 => 'C', 4 => 'G', 8 => 'T', _ => 'N' }
}

fn base_to_acgt_index(base: u8) -> usize {
	match base { b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, _ => 10 }
}

// Function for reading BAM records, with proper user-friendly messages.
// Returns false after reading the last record, or if reading fails.
fn read_bam_record(bam: &mut bam::Reader, record: &mut bam::Record) -> bool {
	match bam.read(record) {
		Ok(true) => true,
		Ok(false) => false,
		Err(bam::Error::TruncatedRecord) => error!("BAM file ended prematurely."),
		Err(bam::Error::InvalidRecord) => error!("Invalid BAM record."),
		Err(e) => error!("Failed reading BAM record: {}", e)
	}
}
