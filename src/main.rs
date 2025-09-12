
#![allow(unused_must_use)]
#![allow(dead_code)]
#![allow(non_snake_case)]

use std::env;
use std::process::exit;

#[macro_use] mod common;
mod call; mod call_targeted; mod analyze; mod counts; mod check_identities;

const USAGE: &str = "
Mutato is a software toolkit for variant calling.

Usage:
  mutato call <genome.fa> <sam/bam/cram>...
  mutato call targeted <genome.fa> <variants.tsv> <sam/bam/cram>...
  mutato analyze <genome.fa> <samples.tsv>
  mutato check identities <mut_files>...
  mutato counts <genome.fa> <sam/bam/cram>...
";

fn main() {
	let args: Vec<String> = env::args().collect();
	let strs: Vec<&str> = args.iter().map(|x| x.as_str()).collect();
	match &strs[..] {
		[_, "call", "targeted", ..] => call_targeted::main(),
		[_, "call", ..] => call::main(),
		[_, "analyze", ..] => analyze::main(),
		[_, "counts", ..] => counts::main(),
		[_, "check", "identities", ..] => check_identities::main(),
		_ => { eprintln!("{}", USAGE); exit(-1); }
	}
}

