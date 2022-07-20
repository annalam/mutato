
use std::env;
use std::process::exit;
use std::io::{Write, stderr};
use docopt::{Docopt, ArgvMap};

#[macro_use] mod common;
mod call;
mod counts;

const USAGE: &str = "
Mutato is a software toolkit for variant calling.

Usage:
  mutato call <genome.fa> <bam_files>...
  mutato counts <genome.fa> <bam_files>...
";

fn main() {
	// TODO: Use match args.as_slice() { [_, "detect"] => ... } after slice
	// pattern matching stabilizes (see issue #23121).
	let args: Vec<String> = env::args().collect();

    if args.len() >= 2 && args[1] == "call" {
    	call::main();
	} else if args.len() >= 2 && args[1] == "counts" {
    	counts::main();
    } else {
    	eprintln!("{}", USAGE); exit(-1);
    }
}

pub fn parse_args(usage: &str) -> ArgvMap {
	Docopt::new(usage).unwrap().parse().unwrap_or_else(|_| {
		writeln!(stderr(), "Invalid arguments.\n{}", usage); exit(-1);
	})
}
