
use std::io::{BufRead, BufReader, stdin};
use std::fs::File;
use std::process::{Command, Stdio};
use std::cmp::Ordering;
use std::collections::HashMap;
use bio::io::fasta;
use docopt::{Docopt, ArgvMap};

macro_rules! error {
	($($arg:tt)+) => ({
		use std::process::exit;
		eprint!("ERROR: "); eprintln!($($arg)+); exit(-1);
	})
}

pub struct FileReader {
	bufread: Box<dyn BufRead>
}

impl FileReader {
	pub fn new(path: &str) -> FileReader {
		let bufread: Box<dyn BufRead> = if path == "-" {
			Box::new(BufReader::new(stdin()))
		} else {
			let file = File::open(path).unwrap_or_else(
				|_| error!("Cannot open file {} for reading.", path));
			if path.ends_with(".gz") {
				Box::new(BufReader::new(Command::new("gunzip").arg("-c")
					.stdout(Stdio::piped()).stdin(file).spawn()
					.unwrap_or_else(|_| error!("Cannot start gunzip process."))
					.stdout.unwrap()))
			} else {
				Box::new(BufReader::new(file))
			}
		};
		FileReader { bufread }
	}

	pub fn read_line(&mut self, line: &mut String) -> bool {
		line.clear();
		match self.bufread.read_line(line) {
			Ok(len) => len > 0,
			_ => { error!("I/O error while reading from file."); }
		}
	}
}

fn chr_num(name: &str) -> u8 {
	let start = if name.starts_with("chr") { 3 } else { 0 };
	if let Ok(num) = name[start..].parse() {
		num
	} else {
		if &name[start..] == "X" { 252 }
		else if &name[start..] == "Y" { 253 }
		else if &name[start..] == "M" { 254 }
		else { 255 }
	} 
}

// Defines a total ordering on chromosome names, in the standard order
pub fn chr_order(a: &str, b: &str) -> Ordering {
	let a_num = chr_num(a);
	let b_num = chr_num(b);
	if a_num == 255 && b_num == 255 {
		// If neither chromosome is a numbered one, sort alphabetically
		a.cmp(&b)
	} else {
		// If one or both chromosomes are numbered, sort by chromosome number.
		a_num.cmp(&b_num)
	}
}


pub struct Genome {
	chrs: Vec<(String, Vec<u8>)>,   // Each chromosome is a name and sequence
	chr_map: HashMap<String, usize>
}

impl Genome {
	pub fn from_fasta(path: &str) -> Genome {
		let reader = fasta::Reader::from_file(&path).unwrap_or_else(|_|
			error!("Could not open genome FASTA file '{}'.", &path));
		let mut chrs: Vec<_> = reader.records().map(|r| {
			let record = r.unwrap();
			let name: String = record.id().into();
			let seq: Vec<u8> = record.seq().iter().map(|x| x.to_ascii_uppercase()).collect();
			(name, seq)
		}).collect();
		chrs.sort_by(|a, b| chr_order(&a.0, &b.0));

		let chr_map: HashMap<_, _> = (0..chrs.len()).map(|c| (chrs[c].0.clone(), c)).collect();
		Genome { chrs, chr_map }
	}

	pub fn name(&self, index: usize) -> &str {
		&self.chrs[index].0
	}

	pub fn sequence(&self, chr: &str) -> Option<&Vec<u8>> {
		match self.chr_map.get(chr) {
			Some(idx) => Some(&self.chrs[*idx].1),
			None => None
		}
	}

	pub fn sequence_by_chr_idx(&self, index: usize) -> &Vec<u8> {
		&self.chrs[index].1
	}

	pub fn chr_idx(&self, chr: &str) -> Option<usize> {
		self.chr_map.get(chr).map(|x| *x)
	}
}

pub fn parse_args(usage: &str) -> ArgvMap {
	Docopt::new(usage).unwrap().parse().unwrap_or_else(|_| {
		error!("Invalid arguments.\n{}", usage);
	})
}
