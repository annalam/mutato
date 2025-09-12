use std::cmp::min;
use std::str;
use std::collections::{VecDeque, HashMap};
use rust_htslib::bam::{self, Read, record::Cigar};
use bio::io::fasta;
use rayon::prelude::*;
use crate::common::{parse_args, read_bam_record, TsvReader};
use std::fs::File;
use std::io::Write;

const USAGE: &str = "
Usage:
  mutato counts [options] <genome.fa> <bam_files>...

Options:
  --min-depth=N             Minimum required depth [default: 10]
  --min-mapq=N              Minimum mapping quality [default: 20]
  --min-baseq=N             Minimum base quality [default: 0]
  --max-allele-fraction=F   Maximum variant allele fraction [default: 0.05]
  --min-end-dist=N          Ignore N bases at both read ends [default: 0]
  --threads=N               Maximum number of threads to use [default: 4]
  --out-dir=PATH            Output path [default: ./]
  --blacklist=FILE          TSV file of blacklisted variants [default: none]

Calculates base substitution error rate across all genomic positions, relative
to a provided reference genome. A blacklist of known mutations and SNPs, as well
as a maximum allowed allele fraction for errors can be provided to mitigate
error rate inflation due to genuine somatic or germline variants.
";


struct SubstitutionMatrix {
    matrix: [[u32; 4]; 4]
}

impl SubstitutionMatrix {
	fn add(&mut self, ref_base: usize, seen_base: usize, n: u32) {
		// ref_base and var_base are given as 0123 = ACGT
		self.matrix[ref_base][seen_base] += n;
	}
    pub fn new() -> SubstitutionMatrix {
		SubstitutionMatrix { matrix: [[0; 4]; 4] }
	}
}

#[derive(Clone)]
struct Pileup {
	total: u32,
	acgt: [u32; 4]
}

type Blacklist = HashMap<String, Vec<u32>>;

#[derive(Default)]
pub struct Settings {
	pub out_dir: String,
	pub min_depth: u32,
	pub min_mapq: u8,
	pub min_baseq: u8,
	pub max_allele_frac: f32,
	pub max_threads: usize,
	pub read_strand: bool,
	pub min_end_distance: u8,
	pub debug: bool
}

pub fn main() {
	let args = parse_args(USAGE);
	let genome_path = args.get_str("<genome.fa>");
	let bam_paths = args.get_vec("<bam_files>");
	let blacklist = args.get_str("--blacklist");

	let mut settings = Settings::default();
	settings.out_dir = args.get_str("--out-dir").into();
	settings.min_depth = args.get_str("--min-depth").parse().unwrap_or_else(
		|_| error!("--min-depth must be a positive integer"));
	settings.min_mapq = args.get_str("--min-mapq").parse().unwrap_or_else(
		|_| error!("--min-mapq must be a positive integer"));
	settings.min_baseq = args.get_str("--min-baseq").parse().unwrap_or_else(
		|_| error!("--min-baseq must be a positive integer"));
	settings.max_allele_frac = args.get_str("--max-allele-fraction").parse().unwrap_or_else(|_| error!("--max-allele-fraction must be a decimal number between 0 and 1"));
	settings.min_end_distance = args.get_str("--min-end-dist").parse()
		.unwrap_or_else(|_| error!("--min-end-dist must be a positive integer."));
	settings.max_threads = args.get_str("--threads").parse().unwrap_or_else(
		|_| error!("--threads must be a positive integer."));

	eprintln!("Reading reference genome into memory...");
	let genome_reader = fasta::Reader::from_file(&genome_path).
		unwrap_or_else(|_| error!("Could not open genome FASTA file '{}'.",
		&genome_path));
	let mut genome: HashMap<String, Vec<u8>> = HashMap::new();
	for r in genome_reader.records() {
		let record = r.unwrap();
		genome.insert(record.id().into(),
			record.seq().iter().map(|x| x.to_ascii_uppercase()).collect());
	}

	let blacklist = if blacklist == "none" { HashMap::new() } else {
		eprintln!("Reading blacklist into memory...");
		read_blacklist(&blacklist)
	};

	// Initialize the thread pool for analyzing multiple BAM files in parallel
	rayon::ThreadPoolBuilder::new().num_threads(settings.max_threads)
		.build_global().unwrap();

	eprintln!("Analyzing BAM files:");
	bam_paths.par_iter().for_each(|bam_path| {
		eprintln!("- {}", &bam_path);
		let bam = bam_path.split("/").last().unwrap();
		let substitution_matrix = analyze_bam(&bam_path, &genome, &blacklist, &settings);
    
		let sample: Vec<_> = bam.split(".").collect();
		let sample = sample[0..sample.len()-1].join("");
		let out = &format!("{}/{}-counts.tsv", &settings.out_dir, sample);
		let mut file = File::create(out).unwrap_or_else(
			|_| error!("Could not create file {}", out));
		
		// Write base substitution matrix
		let ACGT = ['A', 'C', 'G', 'T'];
		for ref_base in 0..4 {
			for seen_base in 0..4 {
				write!(file, "{} -> {}\t{}\n", ACGT[ref_base], ACGT[seen_base],
					substitution_matrix.matrix[ref_base][seen_base]);
			}
		}
	});
}


fn analyze_bam(bam_path: &str, genome: &HashMap<String, Vec<u8>>, 
	blacklist: &Blacklist, settings: &Settings) -> SubstitutionMatrix {

	// Open the BAM file for reading
	let mut bam = bam::Reader::from_path(bam_path).unwrap_or_else(
		|_| error!("Could not open BAM file '{}'.", &bam_path));
	let header = bam.header().clone();
	let chr_names: Vec<&str> = header.target_names().iter().map(
        |x| str::from_utf8(x).unwrap()).collect();
	let empty_list: Vec<u32> = Vec::new();

	let mut pileups: VecDeque<Pileup> = VecDeque::with_capacity(1000);
	let mut curr_chr: u32 = u32::MAX;
	let mut curr_pos: u32 = 0;
	let mut chr_seq = &genome["chr1"];
	let mut chr_black_pos = &empty_list;
	let mut curr_chr_black: usize = 0; // Rolling index to the above chr_black_pos
	let mut substitution_matrix = SubstitutionMatrix::new();

	let mut read = bam::Record::new();
	while read_bam_record(&mut bam, &mut read) {
		if read.is_unmapped() || read.is_duplicate() { continue; }
		if read.is_secondary() || read.is_supplementary() { continue; }
		if read.is_quality_check_failed() { continue; }
		if read.tid() < 0 { error!("Invalid TID < 0."); }
		if read.mapq() < settings.min_mapq { continue; }

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
			
			while curr_chr_black < chr_black_pos.len() &&
				chr_black_pos[curr_chr_black] < curr_pos {
				curr_chr_black += 1;
			}
			
			if (curr_chr_black != chr_black_pos.len() && 
				chr_black_pos[curr_chr_black] == curr_pos) ||
				pileup.total < settings.min_depth {
				curr_pos += 1;
				continue
			}

			let ref_base = base_to_acgt_index(chr_seq[curr_pos as usize - 1].to_ascii_uppercase());
			for base in 0..4 {
				if ref_base == 10 { continue }  // Unknown reference base
				if ref_base != base && pileup.acgt[base] as f32 / pileup.total as f32 > settings.max_allele_frac { continue }
				substitution_matrix.add(ref_base, base, pileup.acgt[base]);
			}

			curr_pos += 1;
		}

		// It is critical that these are updated only *after* reporting
		// the pileups.
		if chr != curr_chr {
			let chr_name = chr_names[chr as usize];
			chr_seq = genome.get(chr_name).unwrap_or_else(||
				error!("BAM file {} refers to region {}, but no such region is found in the provided genome FASTA file.", 
                    &bam_path, &chr_name));
			chr_black_pos = blacklist.get(chr_name).unwrap_or(&empty_list);
			curr_chr_black = 0;
			curr_chr = chr;
		}
		curr_pos = pos;

		// At this point pileups[0] represents chromosome position "curr_pos".
		// We add the read to the pileups vector.
		add_read_to_pileups(&mut pileups, &read, &settings);
	}

	substitution_matrix
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
	read_span += 1;   // We add one extra in case the CIGAR ends in insertion

	// Extend the pileup vector if necessary.
	while pileups.len() < read_span {
		pileups.push_back(Pileup { acgt: [0; 4], total: 0 });
	}

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
					let acgt_idx = match seq.encoded_base(seq_idx) {
						1 => 0,   // A
						2 => 1,   // C
						4 => 2,   // G
						8 => 3,   // T
						_ => 4    // N
					};

					let end_distance = min(seq_idx, seq.len() - 1 - seq_idx) as u32;

					if acgt_idx < 4 && qual[seq_idx] >= settings.min_baseq && end_distance >= settings.min_end_distance as u32 {
						pileups[pileup_idx].total += 1;
						pileups[pileup_idx].acgt[acgt_idx] += 1;
					}
					pileup_idx += 1;
					seq_idx += 1;
				}
			},
			Cigar::Ins(len) => { seq_idx += len as usize; },
			Cigar::Del(len) => { pileup_idx += len as usize; },
			Cigar::RefSkip(len) => { pileup_idx += len as usize; },
			Cigar::SoftClip(len) => { seq_idx += len as usize; },
			Cigar::HardClip(_) => {},
			_ => error!("Unsupported CIGAR element")
		}
	}
}

fn read_blacklist(file: &str) -> Blacklist {
	let mut blacklist = Blacklist::new();
	let file = TsvReader::new(file);
	for cols in file {
		if cols.len() < 2 { 
			error!("Line has < 2 columns in the blacklist file:\n{:?}", cols);
		}
		let positions = blacklist.entry(cols[0].clone()).or_insert_with(|| Vec::new());
		positions.push(cols[1].parse().unwrap_or_else(|_| error!("Invalid position '{}' in blacklist file.\n", cols[1])));
	}
	for pos in blacklist.values_mut() { pos.sort(); pos.dedup(); }
	blacklist
}

fn base_to_acgt_index(base: u8) -> usize {
	match base { b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, _ => 10 }
}

