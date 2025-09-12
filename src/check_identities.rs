use std::str;
use std::collections::HashMap;
use crate::common::{parse_args, TsvReader};

const USAGE: &str = "
Usage:
  mutato check identities [options] <mut_files>...
";

struct Sample {
	name: String,
	patient: String,
	snp_af: Vec<(usize, usize)>
}

pub fn main() {
	let args = parse_args(USAGE);
	let mut_files = args.get_vec("<mut_files>");

	let mut samples: Vec<Sample> = Vec::new();
	let mut snps: HashMap<String, usize> = HashMap::new();

	eprintln!("Reading SNP information into memory...");
	for mut_file in &mut_files {
		let patient = mut_file.replace(".mut", "");
		let mut tsv = TsvReader::new(&mut_file);
		let headers = tsv.read().unwrap();
		let mut patient_samples: Vec<_> = headers[5..headers.len()-1].iter().map(|s|
			Sample { name: s.clone(), patient: patient.clone().into(), snp_af: Vec::new() }).collect();
		for row in tsv {
			if row.len() < 4 + patient_samples.len() { 
				error!("{} contains incomplete line:\n{:?}\n", &mut_file, &row);
			}
			if row[4] != "Germline (ID)" { continue; }
			let snp_id = format!("{}:{}:{}>{}", &row[0], &row[1], &row[2], &row[3]);
			let S = snps.len();
			let r = *snps.entry(snp_id).or_insert(S);
			for s in 0..patient_samples.len() {
				let mut fields = row[5+s].split(':');
				let alt: usize = fields.next().unwrap().parse().unwrap();
				let total: usize = fields.next().unwrap().parse().unwrap();
				patient_samples[s].snp_af.resize(snps.len(), (0, 0));
				patient_samples[s].snp_af[r] = (alt, total);
			}
		}

		samples.extend(patient_samples);
	}

	// Make sure all samples have the same size SNP vector
	for s in &mut samples { s.snp_af.resize(snps.len(), (0, 0)); }

	eprintln!("Comparing SNP profiles...");
	for a in 0..samples.len()-1 {
		for b in a+1..samples.len() {
			let mut concordant = 0;
			let mut informative = 0;
			for r in 0..snps.len() {
				let (a_alt, a_total) = samples[a].snp_af[r];
				let (b_alt, b_total) = samples[b].snp_af[r];
				let a_gt = genotype(a_alt, a_total);
				let b_gt = genotype(b_alt, b_total);
				if a_gt != 0 && b_gt != 0 {
					informative += 1;
					concordant += if a_gt == b_gt { 1 } else { 0 };
				}
			}

			if informative == 0 {
				println!("{}\t{}\t???\tNo informative SNPs", samples[a].name, samples[b].name);
				continue;
			}

			let matched = concordant == informative;
			if (samples[a].patient == samples[b].patient) != matched {
				print!("{}\t{}\t{} / {}\t", samples[a].name, samples[b].name, concordant, informative);
				if matched {
					println!("Unexpected matching SNP profile");
				} else {
					println!("Mismatched SNP profile in samples from patient {}", samples[a].patient);
				}
			}
		}
	}
}

fn genotype(alt: usize, total: usize) -> usize {
	let af = (alt as f32) / (total as f32);
	if total < 20 { 0 }   // Uninformative
	else if af <= 0.02 { 1 }
	else if af >= 0.2 && af <= 0.8 { 2 }
	else if af >= 0.98 { 3 }
	else { 0 } // Uninformative
}