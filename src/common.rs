
use std::io::{BufRead, BufReader, stdin};
use std::fs::File;
use std::process::{Command, Stdio};
use docopt::{Docopt, ArgvMap};

macro_rules! error {
	($($arg:tt)+) => ({
		use std::process::exit;
		eprint!("ERROR: "); eprintln!($($arg)+); exit(-1);
	})
}

pub struct FileReader {
	bufread: Box<BufRead>
}

impl FileReader {
	pub fn new(path: &str) -> FileReader {
		let bufread: Box<BufRead> = if path == "-" {
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

	/*pub fn read_ascii_line(&mut self, line: &mut AsciiString) -> bool {
		line.clear();
		let v = unsafe { std::mem::transmute::<&mut AsciiString, &mut Vec<u8>>(line) };
		self.bufread.read_until(b'\n', v).unwrap() > 0
	}*/
}

pub fn parse_args(usage: &str) -> ArgvMap {
	Docopt::new(usage).unwrap().parse().unwrap_or_else(|_| {
		error!("Invalid arguments.\n{}", usage);
	})
}