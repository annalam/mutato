use std::cmp::min;
use std::str;
//use std::io::Write;
use std::collections::{VecDeque, HashMap};
use rust_htslib::bam::{self, Read};
use rust_htslib::bam::record::Cigar;
use bio::io::fasta;
use itertools::Itertools;
use rayon::prelude::*;
// use smartstring::alias::String;    // Small string optimization
use crate::parse_args;

use std::fs::File;
use std::io::{Write, BufRead, BufReader};

const USAGE: &str = "
Usage:
  mutato counts [options] <genome.fa> <bam_files>...

Options:
  --min-depth=N             Minimum number of overlapping reads in a position [default: 10]
  --min-mapq=N              Minimum mapping quality [default: 20]
  --min-baseq=N             Minimum base quality [default: 0]
  --max-allele-fraction=F   Maximum variant allele fraction [default: 0.05]
  --threads=N               Maximum number of threads to use [default: 4]
  --out-dir=PATH            Output path [default: ./]
  --blacklist=FILE          .tsv file of blacklisted variants [default: none]
  --read-strand=BOOL        Report variants from read's strand insted of original strand (mate 1) [default: false]
  --debug=BOOL              Print variants to stdout [default: false]
";


struct SubstitutionMatrix {
    matrix: [[u32; 4]; 4]
}

impl SubstitutionMatrix {
	fn add(&mut self, ref_base: &char, var_base: &char, n: u32, rev: bool) {
		let mut ref_idx = match ref_base {'A'=>0,'C'=>1,'G'=>2,'T'=>3,_=>return};
		let mut var_idx = match var_base {'A'=>0,'C'=>1,'G'=>2,'T'=>3,_=>return};
		if rev {
			ref_idx = 3 - ref_idx;
			var_idx = 3 - var_idx;
		}
		self.matrix[ref_idx][var_idx] += n;
	}
    pub fn new() -> SubstitutionMatrix {
		SubstitutionMatrix { matrix: [[0; 4]; 4] }
	}
}

#[derive(Clone)]
struct Indel {
	sequence: String,    // Prefixed with '+' for insertion, '-' for deletion
	reads: u32,
	strand_reads: u32
}

#[derive(Clone)]
struct Pileup {
	total: u32,
	acgt: [u32; 4],
	acgt_strand: [u32; 4],
	indels: Vec<Indel>
}

#[derive(Clone)]
struct GPos {
	chr_id: u16,
	pos: u32
}
struct Blacklist {
	chr_ids: HashMap<String, u16>,
	list: Vec<GPos>
}
impl Blacklist {
    pub fn new() -> Blacklist {
		Blacklist { chr_ids: HashMap::new(), list: Vec::new() }
	}
}

pub fn main() {
	let args = parse_args(USAGE);
	let genome_path = args.get_str("<genome.fa>");
	let bam_paths = args.get_vec("<bam_files>");
    let out_dir = args.get_str("--out-dir");
    let blacklist = args.get_str("--blacklist");
	let min_depth: u32 = args.get_str("--min-depth").parse().unwrap_or_else(
		|_| error!("--min-depth must be a positive integer"));
	let min_mapq: u8 = args.get_str("--min-mapq").parse().unwrap_or_else(
		|_| error!("--min-mapq must be a positive integer"));
	let min_baseq: u8 = args.get_str("--min-baseq").parse().unwrap_or_else(
		|_| error!("--min-baseq must be a positive integer"));
	let max_allele_frac: f32 = args.get_str("--max-allele-fraction").parse().unwrap_or_else(
		|_| error!("--max-allele-fraction must be a decimal number between 0 and 1"));
	let max_threads: usize = args.get_str("--threads").parse().unwrap_or_else(
		|_| error!("--threads must be a positive integer."));
	let read_strand: bool = args.get_str("--read-strand").parse().unwrap_or_else(
		|_| error!("--read-strand must be true/false."));
	let debug: bool = args.get_str("--debug").parse().unwrap_or_else(
		|_| error!("--debug must be true/false."));

	if debug && bam_paths.len() > 1 {
		eprintln!("WARNING: Variants are not ordered when running debug to multiple samples simultaneously")
	}

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

	// Read common blacklist
	let blacklist: Blacklist = match blacklist {
		"none" => Blacklist::new(),
		_ => {
			eprintln!("Reading blacklist into memory...");
			let blacklist = read_blacklist(&blacklist);
			eprintln!(" - {} blacklisted genome positions loaded", blacklist.list.len());
			blacklist
		}
	};

	// Initialize the thread pool for analyzing multiple BAM files in parallel
	rayon::ThreadPoolBuilder::new().num_threads(max_threads).build_global()
		.unwrap();

	eprintln!("Analyzing BAM files:");
	bam_paths.par_iter().for_each(|bam_path| {
		if debug {
			println!("CHR\tSTART\tSTOP\tREF\tALT\tREADS (VARIANT/TOTAL)\tVARIANT READS (PRIMARY STRAND/TOTAL)");
		}
		eprintln!("- {}", &bam_path);
		let bam = bam_path.split("/").last().unwrap();
		let (substitution_matrix, ins_counts, del_counts) = analyze_bam(&bam_path, &genome, &blacklist, min_depth, min_mapq, min_baseq, max_allele_frac, read_strand, debug);
    
		let sample: Vec<_> = bam.split(".").collect();
		let sample = sample[0..sample.len()-1].join("");
		let out = &format!("{}/{}-counts.tsv", out_dir, sample);
		let mut file = File::create(out).unwrap_or_else(
			|_| error!("Could not create file {}", out));
		
		// Write base substitution matrix
		write!(file, "Base substitution matrix (row: reference, column: variant)\n");
		write!(file, "\tA\tC\tG\tT\n");
		for (r, row) in substitution_matrix.matrix.iter().enumerate() {
			let ref_base = ['A', 'C', 'G', 'T'][r];
			write!(file, "{}\t{}\n", &ref_base, row.iter().join("\t"));
		}
		write!(file, "\n");
		
		write!(file, "Insertion counts per length:\n");
		write!(file, "{}\n", (1..ins_counts.len()+1).join("\t"));
		let ins_string: String = ins_counts.iter().map( |&x| x.to_string() + "\t").collect();
		write!(file, "{}\n", ins_string);
		
		write!(file, "Deletion counts per length:\n");
		write!(file, "{}\n", (1..(del_counts.len()+1)).join("\t"));
		let del_string: String = del_counts.iter().map( |&x| x.to_string() + "\t").collect();
		write!(file, "{}\n", del_string);
	});

	eprintln!("Analysis complete.");
}


fn analyze_bam(bam_path: &str, genome: &HashMap<String, Vec<u8>>, 
	blacklist: &Blacklist, min_depth: u32, min_mapq: u8, min_baseq: u8, 
	max_allele_frac: f32, read_strand: bool, debug: bool) 
	-> (SubstitutionMatrix, [u32; 20], [u32; 20]) {

	// Open the BAM file for reading
	let mut bam = bam::Reader::from_path(bam_path).unwrap_or_else(
		|_| error!("Could not open BAM file '{}'.", &bam_path));
	let header = bam.header().clone();
	let chr_names: Vec<&str> = header.target_names().iter().map(
        |x| str::from_utf8(x).unwrap()).collect();

	let mut pileups: VecDeque<Pileup> = VecDeque::with_capacity(1000);
	let mut curr_chr: u32 = u32::MAX;
	let mut curr_pos: u32 = 0;
	let mut chr_seq = &genome["chr1"];
	let mut chr_black_pos: Vec<u32> = Vec::new();
	let mut curr_chr_black: usize = 0;
	let mut substitution_matrix = SubstitutionMatrix::new();
	let mut ins_counts: [u32; 20] = [0; 20];
	let mut del_counts: [u32; 20] = [0; 20];

	let mut read = bam::Record::new();
	while read_bam_record(&mut bam, &mut read) {
		if read.is_unmapped() || read.is_duplicate() { continue; }
		if read.is_secondary() || read.is_supplementary() { continue; }
		if read.is_quality_check_failed() { continue; }
		if read.tid() < 0 { error!("Invalid TID < 0."); }
		if read.mapq() < min_mapq { continue; }

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
			
			if curr_chr_black != chr_black_pos.len() && 
				chr_black_pos[curr_chr_black] == curr_pos ||
				pileup.total < min_depth	{
						curr_pos += 1;
				break
			}

			let ref_base = (chr_seq[curr_pos as usize - 1] as char).to_ascii_uppercase();
			const ACGT: [char; 4] = ['A', 'C', 'G', 'T'];
			let mut valid_allele_frac = true;
			for base in 0..4 {
				if ACGT[base] == ref_base { continue; };
				if (pileup.acgt[base] / pileup.total) as f32 > max_allele_frac {
						valid_allele_frac = false;
						break
				}
			}
			if valid_allele_frac {
				for base in 0..4 {
					if pileup.acgt[base] == 0 { continue; };
					if debug && ACGT[base] != ref_base {
						let chr_name = chr_names[chr as usize];
						print!("{}\t{}\t{}\t{}\t{}", &chr_name, curr_pos-1, curr_pos, &ref_base, &ACGT[base]);
						let vaf_percent: f32 = pileup.acgt[base] as f32 / pileup.total as f32 * 100.0;
						print!("\t{}/{} ({:.3}%)", &pileup.acgt[base], &pileup.total, &vaf_percent);
						let strand_percent = pileup.acgt_strand[base] as f32 / pileup.acgt[base] as f32 * 100.0;
						println!("\t{}/{} ({:.3}%)", &pileup.acgt_strand[base], &pileup.acgt[base], &strand_percent);
					}
					substitution_matrix.add(&ref_base, &ACGT[base], pileup.acgt_strand[base], false);
					substitution_matrix.add(&ref_base, &ACGT[base], pileup.acgt[base] - pileup.acgt_strand[base], true);
				}
			}
			
			for indel in &pileup.indels {
				if indel.sequence.len()-1 > ins_counts.len() { continue; };
				let sign = indel.sequence.as_bytes()[0];
				if sign == "+".as_bytes()[0] {
					ins_counts[indel.sequence.len()-2] += indel.reads;
					if debug {
						let chr_name = chr_names[chr as usize];
						let ins_seq = &indel.sequence[1..];
						print!("{}\t{}\t{}\t{}\t{}{}", &chr_name, curr_pos-1, curr_pos, &ref_base, &ref_base, &ins_seq);
						let strand_percent: f32 = indel.strand_reads as f32 / indel.reads as f32 * 100.0;
						println!("\t\t{}/{} ({:.3}%)", &indel.strand_reads, &indel.reads, &strand_percent);
					}
				} else {
					del_counts[indel.sequence.len()-2] += indel.reads;
					if debug {
						let chr_name = chr_names[chr as usize];
						let ref_seq = str::from_utf8(&chr_seq[(curr_pos - 1) as usize .. ((curr_pos - 1) as usize + indel.sequence.len())]).unwrap();
						print!("{}\t{}\t{}\t{}\t{}", &chr_name, curr_pos-1, curr_pos, &ref_seq, &ref_base);
						let strand_percent: f32 = indel.strand_reads as f32 / indel.reads as f32 * 100.0;
						println!("\t\t{}/{} ({:.3}%)", &indel.strand_reads, &indel.reads, &strand_percent);
					}
				}
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
			if blacklist.chr_ids.contains_key(&chr_name.to_string()) {
				let id = blacklist.chr_ids[&chr_name.to_string()];
				chr_black_pos = blacklist.list.iter().filter(|gp| gp.chr_id == id).map(|gp| gp.pos).collect();
			} else {
				chr_black_pos = Vec::new();
			};
			curr_chr_black = 0;
			curr_chr = chr;
		}
		curr_pos = pos;

		// At this point pileups[0] represents chromosome position "curr_pos".
		// We add the read to the pileups vector.
		add_read_to_pileups(&mut pileups, &read, &min_baseq, read_strand);
	}

	(substitution_matrix, ins_counts, del_counts)
}



fn add_read_to_pileups(pileups: &mut VecDeque<Pileup>, read: &bam::Record, min_baseq: &u8, read_strand: bool) {
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
		pileups.push_back(Pileup { acgt: [0; 4], acgt_strand: [0; 4], total: 0, indels: Vec::new() });
	}

	// These variables keep track of our current position
	// within the pileup track and the read sequence, as we
	// read through the CIGAR string.
	let mut pileup_idx = 0;
	let mut seq_idx = 0;
	let mut prior_n = false;

	let seq = read.seq();
	let qual = read.qual();

	// seq represents the strand to be reported (mate 1 or read depending on --read-strand)
	let primary_strand = (read_strand && !read.is_reverse()) || 
							(!read_strand && (!read.is_paired() || read.is_first_in_template()));

	for s in read.cigar().iter() {
		match *s {
			Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
				for _ in 0..len {
					if &qual[seq_idx] >= min_baseq {
						let acgt_idx = match seq.encoded_base(seq_idx) {
							1 => 0,   // A
							2 => 1,   // C
							4 => 2,   // G
							8 => 3,   // T
							_ => 4    // N
						};
						if acgt_idx >= 4 {
							prior_n = true;  // Ambiguous nucleotide
						} else {
							prior_n = false;
							pileups[pileup_idx].total += 1;
							pileups[pileup_idx].acgt[acgt_idx] += 1;
							if primary_strand {
								pileups[pileup_idx].acgt_strand[acgt_idx] += 1;
							};
						};
					}
					pileup_idx += 1;
					seq_idx += 1;
				}
			},
			Cigar::Ins(len) => {
				if pileup_idx == 0 {
					error!("CIGAR strings starting with insertion are not supported.");
				} 

				let mut allele = String::new();
				allele.push('+');     // Mark insertions with a plus
				for _ in 0..len {
					let b = match seq.encoded_base(seq_idx) { 
						1 => 'A', 
						2 => 'C', 
						4 => 'G', 
						8 => 'T', 
						_ => 'N' 
					};
					allele.push(b);
					seq_idx += 1;
				}
				if allele.contains('N') { continue; }  // No ambiguous

				let pileup = &mut pileups[pileup_idx - 1];
				count_indel(pileup, &allele, primary_strand);

				// Indels are assigned to the coordinate of the previous
				// matching base. If the previous matching base was N,
				// it was not added to pileup.total and must be added here.
				if prior_n { pileup.total += 1; }

			},
			Cigar::Del(len) => { 
				if pileup_idx == 0 {
					error!("CIGAR strings starting with deletion are not supported.");
				}

				let mut pileup = &mut pileups[pileup_idx - 1];
				let mut allele = String::new();
				allele.push('-');
				for _ in 0..len {
					allele.push('N');
					pileup_idx += 1;
				}

				count_indel(&mut pileup, &allele, primary_strand);

				// Indels are assigned to the coordinate of the previous
				// MATCHing base. If the previous MATCHing base was an N,
				// it has not been added to pileup.total and must be
				// added here.
				if prior_n { pileup.total += 1; }
			},
            // Soft and hard clips must come last in a CIGAR string.
			// So if we see them, we are done with the read.
			Cigar::SoftClip(_) | Cigar::HardClip(_) => { break; },
			Cigar::RefSkip(len) => { pileup_idx += len as usize; },
			_ => error!("Unsupported CIGAR element")
		}
	}
}


fn read_blacklist(file: &str) -> Blacklist {
	let mut blacklist: Blacklist = Blacklist::new();
	let io = File::open(file).unwrap_or_else(
		|_| error!("Error reading blacklist file: {}", &file));
	let io = BufReader::new(io);
	let mut chr_id_counter: u16 = 0;
	for (line_number, line) in io.lines().enumerate() {
		let line = line.unwrap_or_else(
			|_| error!("Error reading line {} of blacklist file {}", 
						&line_number, &file));
		let line: Vec<_> = line.split("\t").collect();
		if line.len() < 2 { 
			error!("Too few elements in line {} of blacklist file {}", 
					&line_number, &file);
		}
		let chr = line[0];
		if !blacklist.chr_ids.contains_key(&chr.to_string()) {
			blacklist.chr_ids.insert(chr.to_string(), chr_id_counter);
			chr_id_counter += 1;
		}
		blacklist.list.push(GPos {
			chr_id: blacklist.chr_ids[&chr.to_string()],
			pos: line[1].parse().unwrap_or_else(
				|_| error!("Error reading position of variant in line {} of blacklist file {}", 
							&line_number, &file))
		});
	}
	blacklist
}

fn count_indel(pileup: &mut Pileup, seq: &String, primary_strand: bool) {
	for i in 0..pileup.indels.len() {
		let mut indel = &mut pileup.indels[i];
		if &*indel.sequence == seq { // Why derefrence?
			indel.reads += 1;
			if primary_strand {
				indel.strand_reads += 1;
			}
			return;
		}
	}
	let indel = Indel {
		sequence: seq.parse().unwrap(), 
		reads: 1,
		strand_reads: primary_strand as u32
	};
	pileup.indels.push(indel);
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

