use std::cmp::min;
use std::str;
use std::fs::File;
use std::io::{Write, BufWriter};
use std::collections::{HashSet, VecDeque};
use rust_htslib::bam::{self, Read, record::Cigar};
use itertools::Itertools;
use smallvec::SmallVec;
use rayon::prelude::*;
use crate::common::{parse_args, Genome, TsvReader, is_file, is_dir, last_modified, read_bam_record};

const USAGE: &str = "
Usage:
  mutato analyze [options] <genome.fa> <samples.tsv>

Options:
  --alt-reads=N       Minimum alt allele reads to report [default: 5]
  --alt-frac=N        Minimum alt allele fraction to report [default: 0.05]
  --min-mapq=N        Minimum average MAPQ of mutant reads [default: 0]
  --min-end-dist=N    Minimum distance of mutation from read ends [default: 0]
  --bg-ratio=N        Allele fraction must be at least N times higher than in
                      the set of control samples [default: 20]
  --germline-ratio=N  Allele fraction must be at least N times higher than in
                      the matched germline sample [default: 3]

  --min-baseq=N       Ignore bases with a base quality < N [default: 0]
  --ignore-ends=N     Ignore N bases at both read ends [default: 0]
  --max-clip=N        Ignore reads have > N clipped bases [default: 999999]

  --concordant        Only count concordantly aligned read pairs
  --ffpe              Suppress formaldehyde induced artifacts (e.g. FFPE sample)
  --oxidative         Suppress oxidative damage (i.e. oxo-guanine)
  --use-duplicates    Include duplicate reads in analysis

  --hotspots=FILE     Apply relaxed thresholds to mutation hotspots in FILE
  --id-snps=FILE      SNPs to use for identity verification
  --threads=N         Maximum number of threads to use [default: 4]
  --max-patients=N    Maximum patients to analyze in one run [default: 1000]
  --output-dir=DIR    Directory for output files [default: ./mutations]

Variant calling thresholds can be overridden on a per-sample basis by
specifying them in the sample sheet.

A set of control samples can be provided in the sample sheet. These samples
will be used to estimate the background error rate at each genomic
position (separately for each possible mutant allele).

Identified somatic mutations are written into SAMPLE_somatic.tsv files in the
output directory. Germline variants are written into a SAMPLE_germline.tsv
file.
";

// The allele enum is 16 bytes in size (was 32 bytes without the Box).
// Could be 8 bytes with some extra niche optimization.
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum Allele {
	#[default] N, A, C, G, T, Insertion(Box<String>), Deletion(usize)
}
use Allele::*;

#[derive(Default, Clone)]
struct AlleleReads {
	allele: Allele,
	reads: u32,
	sidedness: u32,     // Sum of distances from nearest read end
	mapq: u32,          // Sum of MAPQ values of reads supporting this allele
	read_strand: u32,   // Supporting reads in + strand
	phys_strand: u32    // Supporting reads from ssDNA fragments in + strand
}

#[derive(Default, Clone)]
struct Pileup { alleles: SmallVec<[AlleleReads; 3]> }

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Mutation {
	pub chr: u16,
	pub position: u32,  // Insertion right of splice, or first deleted base
	pub allele: Allele
}

#[derive(Clone, Default)]
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
	pub germline_ratio: f32,
	pub bg_ratio: f32,
	pub min_mapq: u8,
	pub min_end_distance: u8,
	pub max_clip: usize,
	pub min_baseq: u8,
	pub ignore_ends: u8,
	pub max_patients: usize,
	pub concordant: bool,
	pub use_duplicates: bool,
	pub suppress_ffpe_artifacts: bool,
	pub suppress_oxidative_damage: bool,
	pub output_dir: String,

	pub hotspots: HashSet<Mutation>,
	pub id_snps: HashSet<Mutation>
}

struct Sample {
	name: String,
	patient: String,
	germline: bool,   // Can sample be used as matched germline control?
	background: bool, // Can sample be included in background error estimation?

	// Sample-specific settings for mutation calling
	alt_frac: f32
}

pub fn main() {
	let args = parse_args(USAGE);
	let max_threads: usize = args.get_str("--threads").parse().unwrap_or_else(
		|_| error!("--threads must be a positive integer."));

	let mut settings = Settings::default();
	settings.alt_reads = args.get_str("--alt-reads").parse().unwrap_or_else(
		|_| error!("--alt-reads must be a positive integer"));
	let alt_frac_global = args.get_str("--alt-frac").parse().unwrap_or_else(
		|_| error!("--alt-frac must be a fraction between 0 and 1."));
	settings.germline_ratio = args.get_str("--germline-ratio").parse().unwrap_or_else(|_| error!("--germline-ratio must be a non-negative number"));
	settings.bg_ratio = args.get_str("--bg-ratio").parse().unwrap_or_else(|_|
		error!("--bg-ratio must be a non-negative number"));
	settings.min_mapq = args.get_str("--min-mapq").parse().unwrap_or_else(|_| error!("--min-mapq must be a non-negative integer."));
	settings.min_end_distance = args.get_str("--min-end-dist").parse().unwrap_or_else(|_| error!("--min-end-dist must be a non-negative number."));
	settings.max_clip = args.get_str("--max-clip").parse().
		unwrap_or_else(|_| error!("--max-clip must be a non-negative number."));
	settings.min_baseq = args.get_str("--min-baseq").parse().unwrap_or_else(
		|_| error!("--min-baseq must be a non-negative integer."));
	settings.ignore_ends = args.get_str("--ignore-ends").parse()
		.unwrap_or_else(|_| error!("--ignore-ends must be a non-negative integer."));
	settings.max_patients = args.get_str("--max-patients").parse().unwrap_or_else(|_| error!("--max-patients must be a positive integer"));
	settings.concordant = args.get_bool("--concordant");
	settings.use_duplicates = args.get_bool("--use-duplicates");
	settings.suppress_ffpe_artifacts = args.get_bool("--ffpe");
	settings.suppress_oxidative_damage = args.get_bool("--oxidative");
	settings.output_dir = args.get_str("--output-dir").into();

	// Check that the output directory exists (so there's no surprise later)
	if is_dir(&settings.output_dir) == false {
		error!("Output directory '{}' does not exist.", &settings.output_dir);
	}

	// Read sample sheet into memory
	let mut sample_sheet = TsvReader::new(args.get_str("<samples.tsv>"));
	let headers = sample_sheet.read().unwrap();
	let patient_col = headers.iter().position(|x| x == "PATIENT").unwrap_or(99);
	let sample_col = headers.iter().position(|x| x == "SAMPLE").unwrap_or(99);
	let type_col = headers.iter().position(|x| x == "TYPE").unwrap_or(99);
	if patient_col == 99 || sample_col == 99 || type_col == 99 {
		error!("Required headers PATIENT, SAMPLE and TYPE were not all found on the first row of the sample sheet.");
	}
	let alt_frac_col = headers.iter().position(|x| x == "ALT_FRAC");
	let min_columns = *[patient_col, sample_col, type_col].iter().max().unwrap();
	let mut samples: Vec<Sample> = Vec::new();
	for row in sample_sheet {
		if row.len() <= min_columns {
			error!("Sample sheet line has too few columns:\n{:?}", &row);
		}
		let patient = &row[patient_col];
		let sample = &row[sample_col];
		let category = row[type_col].to_lowercase();

		// Read sample-specific mutation calling settings, or use the global
		let alt_frac: f32 = match alt_frac_col {
			Some(c) if c < row.len() => row[c].parse().unwrap(),
			_ => alt_frac_global
		};

		// Check that a corresponding BAM file exists in the current folder
		if is_file(&format!("{}.bam", sample)) == false {
			error!("No BAM file found for sample {}.", &sample);
		}

		samples.push(Sample { name: sample.clone(), patient: patient.clone(),
			germline: category.contains("wbc") || category.contains("gdna") || category.contains("germline"),
			background: category.contains("background"),
			alt_frac });
	}

	if samples.iter().all(|s| s.background == false) {
		error!("The cohort does not contain any background control samples.");
	}

	eprintln!("Sample sheet contains {} samples ({} germline, {} background controls).", samples.len(),
		samples.iter().filter(|s| s.germline).count(),
		samples.iter().filter(|s| s.background).count());

	// Determine which patients in the sample sheet require a new analysis.
	// We only redo an analysis if the set of samples for a patient has changed.
	// We don't redo an analysis if the set of control samples has changed.
	let mut patients: Vec<String> = samples.iter().map(|s| s.patient.clone()).collect();
	patients.sort();
	patients.dedup();

	let mut patients_to_analyze: Vec<String> = Vec::new();
	for patient in &patients {
		let newest_bam_mtime = samples.iter().filter(|s| s.patient == *patient)
			.map(|s| last_modified(&format!("{}.bam", &s.name)))
			.max().unwrap();

		let mut_file = format!("{}/{}.mut", &settings.output_dir, &patient);
		if !is_file(&mut_file) || last_modified(&mut_file) < newest_bam_mtime {
			patients_to_analyze.push(patient.clone());
		}
		// TODO: Also check that analysis thresholds haven't changed?
	}
	if patients_to_analyze.is_empty() {
		eprintln!("All patient *.mut files are up-to-date.");
		return;
	}

	if patients.len() > patients_to_analyze.len() {
		eprintln!("{} / {} patients are up-to-date and do not require a re-analysis.", patients.len() - patients_to_analyze.len(),
			patients.len());
	}

	if patients_to_analyze.len() > settings.max_patients {
		eprintln!("Restricting analysis to {} patients to prevent excessive memory use (see --max-patients=N):", settings.max_patients);
		// Prioritize patients that have cancer samples
		patients_to_analyze.sort_by_key(|p| if samples.iter().any(|s| s.patient == *p && s.germline == false) { 1 } else { 2 });
		patients_to_analyze.truncate(settings.max_patients);
		for p in &patients_to_analyze { eprintln!("- {}", p); }
	}

	// Build a list of samples for the patients that require a re-analysis.
	// Since the goal of the first analysis pass is to identify candidate
	// variants, we don't need to include samples that are only used
	// for calculating the background error rate.
	let first_pass_samples: Vec<&Sample> = samples.iter().filter(|s| patients_to_analyze.contains(&s.patient)).collect();

	eprintln!("Reading reference genome into memory...");
	let genome = Genome::from_fasta(&args.get_str("<genome.fa>"));

	// Read hotspot mutation into memory
	let hotspots_path = args.get_str("--hotspots");
	if hotspots_path.is_empty() == false {
		eprintln!("Reading list of hotspot mutations...");
		let tsv = TsvReader::new(&hotspots_path);
		for row in tsv {
			if row.len() < 3 { error!("Line has too few columns:\n{:?}", &row);}
			if row[2].len() != 1 { unimplemented!(); }
			settings.hotspots.insert(Mutation {
				chr: genome.chr_idx(&row[0]) as u16,
				position: row[1].parse().unwrap(),
				allele: ascii_to_allele(row[2].as_bytes()[0])
			});
		}
	}

	// Read list of patient identity verification SNPs into memory
	let id_snps_path = args.get_str("--id-snps");
	if id_snps_path.is_empty() == false {
		eprintln!("Reading list of patient identity verification SNPs...");
		let tsv = TsvReader::new(&id_snps_path);
		for row in tsv {
			if row.len() < 3 { error!("Line has too few columns:\n{:?}", &row);}
			if row[2].len() != 1 { unimplemented!(); }
			settings.id_snps.insert(Mutation {
				chr: genome.chr_idx(&row[0]) as u16,
				position: row[1].parse().unwrap(),
				allele: ascii_to_allele(row[2].as_bytes()[0])
			});
		}
	}

	// Initialize the thread pool for analyzing multiple BAM files in parallel
	rayon::ThreadPoolBuilder::new().num_threads(max_threads).build_global()
		.unwrap();

	eprintln!("Analyzing BAM files for candidate variants:");
	let mut sample_variants: Vec<_> = first_pass_samples.par_iter().map(|&sample| {
		eprintln!("- {}", &sample.name);
		let bam_path = format!("{}.bam", &sample.name);
		let mut variants = analyze_bam_firstpass(&bam_path, &genome, &settings, &sample);
		variants.sort();  // Sort in preparation for k-way merge
		variants
	}).collect();

	// We take the union of candidate variants across all samples, and include
	// identity verification SNPs if specified by the user.
	let mut id_snps_sorted: Vec<_> = settings.id_snps.iter().cloned().collect();
	id_snps_sorted.sort();
	sample_variants.push(id_snps_sorted);
	let variants: Vec<Mutation> = itertools::kmerge(sample_variants).dedup().collect();

	eprintln!("Found {} candidate variants.", variants.len());

	// In the second pass, we will analyze all samples from the patients
	// selected for this analysis, as well as all background control samples
	// in the entire cohort.
	let second_pass_samples: Vec<&Sample> = samples.iter().filter(|s| s.background || patients_to_analyze.contains(&s.patient)).collect();

	eprintln!("Calculating evidence matrix...");
	let evidence: Vec<_> = second_pass_samples.par_iter().map(|&sample| {
		let bam_path = format!("{}.bam", &sample.name);
		analyze_bam_secondpass(&bam_path, &genome, &variants, &settings)
	}).collect();

	eprintln!("Writing mutation calls to disk...");
	patients_to_analyze.par_iter().for_each(|patient| {
		report_variants(&patient, &second_pass_samples, &variants, &evidence, &settings, &genome);
	});
}


fn analyze_bam_firstpass(bam_path: &str, genome: &Genome, settings: &Settings, sample: &Sample) -> Vec<Mutation> {

	// Apply a more stringent allele fraction threshold for germline controls
	let min_maf = if sample.germline { 0.1 } else { sample.alt_frac };
	let hotspot_maf = min_maf / 2.0;

	// For hotspot mutations we halve the required number of mutant allele reads
	// (rounded up, hence the +1). TODO: Use div_ceil() once stabilized.
	let hotspot_alt_reads = (settings.alt_reads + 1) / 2;

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
		if read.is_unmapped() { continue; }
		if settings.use_duplicates == false && read.is_duplicate() { continue; }
		if read.is_secondary() || read.is_supplementary() { continue; }
		if read.is_quality_check_failed() { continue; }
		if settings.concordant && read.is_paired() {
			if read.tid() != read.mtid() { continue; }
		}

		let chr = read.tid() as u32;
		let pos = read.pos() as u32 + 1;   // Convert to one-based coordinate

		// Check how many pileups we have finalized (i.e. will no longer change)
		let to_report = if chr != curr_chr {
			pileups.len()
		} else if pos > curr_pos {
			min((pos - curr_pos) as usize, pileups.len())
		} else if pos < curr_pos {
			error!("{} is not sorted by position.", &bam_path)
		} else {
			0
		};

		// Go through the finalized pileups, and check if any meet our
		// thresholds.
		for pileup in pileups.drain(0..to_report) {
			let ref_allele = match chr_seq.get(curr_pos as usize - 1) {
				Some(&b) => ascii_to_allele(b), None => N
			};

			for a in &pileup.alleles {
				if a.allele == ref_allele || a.allele == N { continue; }
				if a.reads < hotspot_alt_reads { continue; }

				let maf = (a.reads as f32) / (pileup.total_reads(&a.allele) as f32);
				if maf < hotspot_maf { continue; }

				// Now that we know we reach at least the hotspot thresholds,
				// we check if this is a known mutation hotspot. If it is not,
				// we apply the non-relaxed thresholds.
				let chr_idx = genome.chr_idx(&chr_names[curr_chr as usize]);
				let mutation = Mutation { chr: chr_idx as u16, position: curr_pos, allele: a.allele.clone() };
				let hotspot = settings.hotspots.contains(&mutation);
				if hotspot == false {
					if a.reads < settings.alt_reads { continue; }
					if maf < min_maf { continue; }
				}
				mutations.push(mutation);
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
	let mut chr_muts: VecDeque<usize> = VecDeque::new();
	let mut evidence = vec![Evidence { alt_reads: 0, total_reads: 0, mapq: 0, sidedness: 0, strand: 0, other_strand: 0 }; mutations.len()];

	let mut read = bam::Record::new();
	while read_bam_record(&mut bam, &mut read) {
		if read.is_unmapped() { continue; }
		if settings.use_duplicates == false && read.is_duplicate() { continue; }
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

				// This genomic position is a candidate variant, so we store 
				// the number of supporting reads in the matrix.
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



fn report_variants(patient: &str, samples: &Vec<&Sample>, variants: &Vec<Mutation>, evidence: &Vec<Vec<Evidence>>, settings: &Settings, genome: &Genome) {

	let out_path = format!("{}/{}.mut", &settings.output_dir, &patient);
	let file = File::create(&out_path).unwrap_or_else(
		|_| error!("Cannot open file {} for writing.", &out_path));
	let mut out = BufWriter::new(file);

	// Determine the set of tumor samples for this patient
	let tumor_samples: Vec<usize> = (0..samples.len()).filter(|&s| samples[s].germline == false && samples[s].patient == patient).collect();

	// Determine the set of matched germline control samples for this patient
	let germline_controls: Vec<usize> = (0..samples.len()).filter(|&s| samples[s].germline && samples[s].patient == patient).collect();

	// Determine the set of background control samples
	let bg_controls: Vec<usize> = (0..samples.len()).filter(|&s| samples[s].background).collect();

	write!(out, "CHROM\tPOSITION\tREF\tALT\tNOTES");
	for s in &tumor_samples { write!(out, "\t{}", &samples[*s].name); }
	for s in &germline_controls { write!(out, "\t{}", &samples[*s].name); }
	write!(out, "\tBACKGROUND%\n");

	for v in 0..variants.len() {

		// Calculate background and germline allele fractions for this variant
		let bg_frac = bg_controls.iter().map(|&s| evidence[s][v].alt_reads as f32).sum::<f32>() / bg_controls.iter().map(|&s| evidence[s][v].total_reads as f32).sum::<f32>();
		let germline_total_reads = germline_controls.iter().map(|&s| evidence[s][v].total_reads).sum::<u32>();
		let germline_frac = germline_controls.iter().map(|&s| evidence[s][v].alt_reads as f32).sum::<f32>() / (germline_total_reads as f32);

		// Check if we are looking at a clear germline variant
		let mut origin: &str = "";  // Somatic, germline, or unknown
		if settings.id_snps.contains(&variants[v]) {
			origin = "Germline (ID)";
		} else if germline_frac >= 0.3 && germline_total_reads >= 10 {
			if germline_frac <= 0.7 {
				origin = "Germline (heterozygous)";
			} else if germline_frac >= 0.95 {
				origin = "Germline (homozygous)";
			}
		}

		// If we are not looking at a clear germline variant, see
		// if we have a somatic mutation or other variant that has significant
		// support relative to the background error rate.
		let mut somatic: Vec<usize> = Vec::new();
		if origin == "" {
			let is_hotspot = settings.hotspots.contains(&variants[v]);

			// Check indel length. We don't apply MAPQ threshold to long indels.
			let indel_len = match &variants[v].allele { Insertion(ins) => ins.len(), Deletion(len) => *len, _ => 0 };

			// Generate a list of tumor samples where the mutation is found
			for &s in &tumor_samples {			
				let e = &evidence[s][v];
				let maf = e.alt_reads as f32 / e.total_reads as f32;

				// Calculate mutation calling thresholds. For hotspot mutations
				// the thresholds are halved (rounded up).
				let min_reads = if is_hotspot { (settings.alt_reads + 1) / 2 } else { settings.alt_reads };
				let min_maf = if is_hotspot { samples[s].alt_frac / 2.0 } else { samples[s].alt_frac }; 

				if e.alt_reads < min_reads { continue }
				if maf < min_maf { continue }

				if is_hotspot == false {
					if e.mapq < settings.min_mapq && indel_len <= 5 { continue }
					if e.sidedness < settings.min_end_distance { continue }
				}

				// If a variant meets all criteria for a somatic mutation but
				// does not exceed the background rate, we still report it out
				// as an "Background" variant as it could be a germline variant
				let min_maf_bg = if is_hotspot { settings.bg_ratio * bg_frac / 2.0 } else { settings.bg_ratio * bg_frac };
				if maf < min_maf_bg {
					origin = "Background";
					continue;
				}

				// We only mark a mutation as "Somatic" if we have a matched
				// germline sample available and the sample has an allele 
				// fraction significantly higher than the germline sample.
				if germline_total_reads >= 5 && maf >= settings.germline_ratio * germline_frac {
					origin = "Somatic";
					somatic.push(s);
				} else if origin == "" {
					origin = "Unknown";
				}
			}
		}

		if origin == "" { continue }

		let chr = variants[v].chr;
		let pos = variants[v].position as usize;
		let chr_seq = genome.sequence_by_chr_idx(chr as usize);
		write!(out, "{}\t", &genome.name(chr as usize));

		match &variants[v].allele {
			// When reporting deletions, we need to subtract 1 from the position
			Insertion(ins) => write!(out, "{}\t{}\t{}{}", pos, chr_seq[pos - 1] as char, chr_seq[pos - 1] as char, ins),
			Deletion(len) => write!(out, "{}\t{}\t{}", pos - 1, str::from_utf8(&chr_seq[pos-2..pos+len-1]).unwrap(), chr_seq[pos - 2] as char),
			_ => write!(out, "{}\t{}\t{}", pos, chr_seq[pos - 1] as char, &variants[v].allele.to_char())
		};
		write!(out, "\t{}", &origin);

		for s in tumor_samples.iter().chain(germline_controls.iter()) {
			let e = &evidence[*s][v];
			write!(out, "\t{}:{}", e.alt_reads, e.total_reads);
			if e.alt_reads == 0 {
				write!(out, ":0::0");
			} else {
				write!(out, ":{}:{},{}:{}", e.mapq, e.strand, e.other_strand, e.sidedness);
				if somatic.contains(&s) { write!(out, ":*"); }
			}
		}
		write!(out, "\t{:.3}\n", bg_frac * 100.0);  // Background
	}
}



fn add_read_to_pileups(pileups: &mut VecDeque<Pileup>, read: &bam::Record, settings: &Settings) {
	// Calculate how long the pileup vector needs to be to accommodate
	// the information from this read.
	let mut read_span = 0;
	let mut total_clip = 0;
	for s in read.cigar().iter() {
		match *s {
			Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) |
			Cigar::Del(len) | Cigar::RefSkip(len) => {
				read_span += len as usize;
			},
			Cigar::Ins(_) => {},
			Cigar::SoftClip(len) | Cigar::HardClip(len) => {
				total_clip += len as usize;
			},
			_ => error!("Unsupported CIGAR element found.")
		}
	}
	read_span += 1;   // We add one extra in case the CIGAR ends in insertion

	if total_clip > settings.max_clip { return }

	// Extend the pileup vector if necessary.
	while pileups.len() < read_span { pileups.push_back(Pileup::default()); }

	// These variables keep track of our current position
	// within the pileup track and the read sequence, as we
	// read through the CIGAR string.
	let mut pileup_idx = 0;
	let mut seq_idx = 0;

	let seq = read.seq();
	let qual = read.qual();

	for s in read.cigar().iter() {
		match *s {
			Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
				for _ in 0..len {
					let allele = bam_base_to_allele(seq.encoded_base(seq_idx));
					if qual[seq_idx] >= settings.min_baseq {
						let pileup = &mut pileups[pileup_idx];
						let end_distance = min(seq_idx, seq.len() - 1 - seq_idx) as u32;
						if end_distance >= settings.ignore_ends as u32 {
							pileup.count_allele(allele, &read, end_distance, read.mapq());
						}
					}
					pileup_idx += 1;
					seq_idx += 1;
				}
			},
			Cigar::Ins(len) => {
				let mut ins = String::new();
				for _ in 0..len {
					ins.push(bam_base_to_char(seq.encoded_base(seq_idx)));
					seq_idx += 1;
				}
				if ins.contains('N') { continue; }  // No ambiguous				

				// Distance from nearest read end
				let end_distance = min(seq_idx - (len as usize),
					seq.len() - 1 - seq_idx) as u32;
				if end_distance < settings.ignore_ends as u32 { continue; }

				// Insertion allele is stored in the pileup structure that
				// is located immediately left of the splice.
				// TODO: Allow a read to begin with an insertion
				// (need to support pileups[-1] indexing here to achieve that).
				if pileup_idx == 0 { continue; }
				let pileup = &mut pileups[pileup_idx - 1];
				let allele = Allele::Insertion(Box::new(ins));
				pileup.count_allele(allele, &read, end_distance, read.mapq());
			},
			Cigar::Del(len) => {
				let allele = Allele::Deletion(len as usize);

				let end_distance = min(seq_idx, seq.len() - 1 - seq_idx) as u32;
				if end_distance < settings.ignore_ends as u32 { continue; }

				// Deletion allele is stored in the pileup structure located
				// atop the first deleted base.
				let pileup = &mut pileups[pileup_idx];
				pileup.count_allele(allele, &read, end_distance, read.mapq());

				pileup_idx += len as usize;
			},
			Cigar::RefSkip(len) => { pileup_idx += len as usize; },
			Cigar::SoftClip(len) => { seq_idx += len as usize; },
			Cigar::HardClip(_) => {},
			_ => error!("Unsupported CIGAR element")
		}
	}
}


fn calculate_mutation_evidence(mutation: &Mutation, pileup: &Pileup,
	genome: &Genome, fix_ffpe: bool, fix_oxog: bool) -> Evidence {
	
	let total_reads = pileup.total_reads(&mutation.allele);
	let total_plus = pileup.total_reads_plus_strand(&mutation.allele);

	// If there were no reads for the mutant allele, generate a zeroed
	// AlleleReads object so the FFPE and oxidative damage calculations
	// can be carried out.
	let dummy = AlleleReads::default();
	let a = pileup.find_allele(&mutation.allele).unwrap_or(&dummy);
	let mut evidence = Evidence {
		alt_reads: a.reads, total_reads,
		mapq: (a.mapq as f32 / a.reads as f32).round() as u8,
		sidedness: (a.sidedness as f32 / a.reads as f32).round() as u8,
		strand: (a.read_strand as f32 / a.reads as f32 * 100.0).round() as u8,
		other_strand: ((total_plus - a.read_strand) as f32 / (total_reads - a.reads) as f32 * 100.0).round() as u8
	};

	// Suppress formaldehyde-induced DNA damage (C>T) and
	// oxidative DNA damage (G>T) by only counting sequenced ssDNA
	// fragments that originate from the strand that cannot be affected
	// by these forms of single-stranded damage.
	if fix_ffpe || fix_oxog {
		let ref_base = genome.sequence_by_chr_idx(mutation.chr as usize)[mutation.position as usize - 1];
		let total_phys_plus = pileup.total_reads_phys_plus_strand(&mutation.allele);

		if fix_ffpe && ref_base == b'C' && mutation.allele == T {
			// Only count ssDNA fragments representing minus strand
			evidence.alt_reads = a.reads - a.phys_strand;
			evidence.total_reads = total_reads - total_phys_plus;
		} else if fix_ffpe && ref_base == b'G' && mutation.allele == A {
			// Only count ssDNA fragments representing plus strand
			evidence.alt_reads = a.phys_strand;
			evidence.total_reads = total_phys_plus;
		} else if fix_oxog && ref_base == b'G' && mutation.allele == T {
			// Only count ssDNA fragments representing minus strand
			evidence.alt_reads = a.reads - a.phys_strand;
			evidence.total_reads = total_reads - total_phys_plus;
		} else if fix_oxog && ref_base == b'C' && mutation.allele == A {
			// Only count ssDNA fragments representing plus strand
			evidence.alt_reads = a.phys_strand;
			evidence.total_reads = total_phys_plus;
		}
	}

	evidence
}

fn bam_base_to_char(encoded: u8) -> char {
	match encoded { 1 => 'A', 2 => 'C', 4 => 'G', 8 => 'T', _ => 'N' }
}

fn bam_base_to_allele(encoded: u8) -> Allele {
	match encoded { 1 => A, 2 => C, 4 => G, 8 => T, _ => N }
}

fn ascii_to_allele(base: u8) -> Allele {
	match base { b'A' => A, b'C' => C, b'G' => G, b'T' => T, _ => N }
}

impl Pileup {

	fn total_reads(&self, mutation: &Allele) -> u32 {

		// The total depth at a given genomic position is calculated differently
		// for base substitutions and indels. For base substitutions, ambiguous
		// N bases are not counted towards the total. For indels, ambiguous N 
		// bases are counted towards the total.
		let mut total = 0;
		for a in &self.alleles {
			// Insertions are not counted because they are actually located
			// between reference genome bases and insertion-containing reads
			// would otherwise be double-counted (since both the
			// insertion allele and its left-side base are recorded here).
			total += match a.allele {
				Insertion(_) => 0,
				N => if mutation.is_indel() { a.reads } else { 0 },
				_ => a.reads
			};
		}
		return total;
	}

	fn total_reads_plus_strand(&self, mutation: &Allele) -> u32 {
		let mut total = 0;
		for a in &self.alleles {
			total += match a.allele {
				Insertion(_) => 0,
				N => if mutation.is_indel() { a.read_strand } else { 0 },
				_ => a.read_strand
			};
		}
		return total;
	}

	fn total_reads_phys_plus_strand(&self, mutation: &Allele) -> u32 {
		let mut total = 0;
		for a in &self.alleles {
			total += match a.allele {
				Insertion(_) => 0,
				N => if mutation.is_indel() { a.phys_strand } else { 0 },
				_ => a.phys_strand
			};
		}
		return total;
	}

	fn find_allele(&self, allele: &Allele) -> Option<&AlleleReads> {
		self.alleles.iter().find(|a| a.allele == *allele)
	}

	fn count_allele(&mut self, allele: Allele, read: &bam::Record, end_distance: u32, mapq: u8) {

		let allele_idx = self.alleles.iter().position(|a| a.allele == allele);
		let a = match allele_idx {
			Some(idx) => &mut self.alleles[idx],
			None => {
				self.alleles.push(AlleleReads { allele, ..Default::default() });
				self.alleles.last_mut().unwrap()
			}
		};

		a.reads += 1;
		a.read_strand += if read.is_reverse() { 0 } else { 1 };

		// In typical Illumina paired end sequencing + library preparation,
		// the mate #1 strand represents the original physical genome strand of
		// the single-stranded DNA fragment that was sequenced. Here we figure
		// out that physical strand.
		let phys_minus = read.is_reverse() == read.is_first_in_template();
		a.phys_strand += if phys_minus { 0 } else { 1 };

		// Running totals, division later
		a.sidedness += end_distance;
		a.mapq += mapq as u32;
	}
}

impl Allele {
	fn is_indel(&self) -> bool {
		match self { Insertion(_) | Deletion(_) => true, _ => false }
	}
	fn to_char(&self) -> char {
		match self { A => 'A', C => 'C', G => 'G', T => 'T', _ => '?' }
	}
}

