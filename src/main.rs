
use std::env;
use std::process::exit;

#[macro_use] mod common;
mod call;
mod call_targeted;
mod counts;

const USAGE: &str = "
Mutato is a software toolkit for variant calling.

Usage:
  mutato call <genome.fa> <bam_files>...
  mutato call targeted <genome.fa> <variants.tsv> <bam_files>...
  mutato counts <genome.fa> <bam_files>...
";

fn main() {
	let args: Vec<String> = env::args().collect();

	// Convert to Vec<&str> to allow slice pattern matching
	let strs: Vec<&str> = args.iter().map(|x| x.as_str()).collect();
	match &strs[..] {
		[_, "call", "targeted", ..] => call_targeted::main(),
		[_, "call", ..] => call::main(),
		[_, "counts", ..] => counts::main(),
		_ => { eprintln!("{}", USAGE); exit(-1); }
	}
}

