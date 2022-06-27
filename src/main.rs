extern crate getopts;
use getopts::Options;
use std::env::args;
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::iter::zip;
use std::path::Path;
use std::process;

fn main() {
	let args: Vec<_> = args().collect();
	let program = args[0].clone();

	let mut opts = Options::new();
	opts.optflag("h", "help", "print this help menu");
	opts.optopt(
		"o",
		"",
		"generate gain, loss and diff files and output them with the prefix NAME",
		"NAME",
	);

	opts.optopt("d",
				"",
				"a location in the tested file is considered to be similar to the one in the target file if the distance between both locations is less than FLOAT fraction of the read length [1.0]",
				"FLOAT",
	);

	opts.optopt(
		"q",
		"",
		"output the results using UINT1,UINT2,... (such as UINT1 > UINT2 > ...) as the quality thresholds [60,10,1,0]",
		"UINT1,UINT2,...",
	);

	opts.optopt(
		"m",
		"",
		"Comparison mode [all]\n \
		- all: \n \
		- prim_ta: \n \
		- prim:",
		"STR",
	);

	let matches = match opts.parse(&args[1..]) {
		Ok(m) => m,
		Err(f) => {
			eprintln!("[ERROR]: {}", f.to_string());
			process::exit(1);
		}
	};

	if matches.opt_present("h") {
		print_usage(&program, opts);
		return;
	}

	let output = matches.opt_str("o");

	let distance: f32 = matches.opt_get_default("d", 1.0).unwrap();

	let qualities: Vec<u8> = match matches.opt_str("q") {
		Some(q) => q.split(',').map(|x| x.parse().unwrap()).collect(),
		None => [60, 10, 1, 0].to_vec(),
	};

	let p_only = matches.opt_present("m");

	if matches.free.len() != 2 {
		print_usage(&program, opts);
		return;
	}

	let path = Path::new(&matches.free[0]);
	let sam_tgt = match parse_sam(&path) {
		Ok(sam_tgt) => sam_tgt,
		Err(err) => {
			eprintln!("[ERROR]: {},", err);
			process::exit(1);
		}
	};

	let path = Path::new(&matches.free[1]);
	let sam_test = match parse_sam(&path) {
		Ok(sam_test) => sam_test,
		Err(err) => {
			eprintln!("[ERROR]: {},", err);
			process::exit(1);
		}
	};

	//sam_tgt.iter().for_each(|x| println!("{}\t", x));
	if sam_tgt.len() != sam_test.len() {
		eprintln!(
			"[ERROR]: Not the same number of reads in each SAM file (target: {} and test: {})",
			sam_tgt.len(),
			sam_test.len()
		);
		process::exit(1);
	}

	eprintln!("[INFO]: {} reads", sam_tgt.len());

	match compare_sam(&sam_tgt, &sam_test, distance, &qualities, &output, p_only) {
		Err(err) => {
			eprintln!("[ERORR]: {}", err);
			process::exit(1);
		}
		_ => {}
	}
}

fn print_usage(program: &str, opts: Options) {
	let brief = format!("Usage: {} [options] <target.sam> <test.sam>", program);
	print!("{}", opts.usage(&brief));
}

#[derive(Clone)]
struct Secondary {
	rname: String,
	//TODO: check if there are also supplementatry alignments for secondaries
	pos: u32,
	strand: bool,
	als: i32,
}

#[derive(Clone)]
struct Mapped {
	qname: String,
	len: u32,
	rname: String,
	pos_min: u32,
	pos_max: u32,
	strand: bool,
	mapq: u8,
	als: i32,
	secondaries: Vec<Secondary>,
}

#[derive(Clone)]
struct Unmapped {
	qname: String,
	len: u32,
}

#[derive(Clone)]
enum Sam {
	Mapped(Mapped),
	Unmapped(Unmapped),
}

impl fmt::Display for Mapped {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		write!(
			f,
			"{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
			self.qname,
			self.len,
			self.rname,
			self.pos_min,
			self.pos_max,
			if self.strand { '-' } else { '+' },
			self.mapq,
			self.als
		)?;
		self.secondaries
			.iter()
			.try_for_each(|x| write!(f, "\t[{}\t{}\t{}\t{}]", x.rname, x.pos, x.strand, x.als))?;
		fmt::Result::Ok(())
	}
}

impl fmt::Display for Unmapped {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		write!(f, "{}\t{}", self.qname, self.len)
	}
}

impl fmt::Display for Sam {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		match self {
			Sam::Mapped(x) => write!(f, "{}", x),
			Sam::Unmapped(x) => write!(f, "{}", x),
		}
	}
}

fn parse_sam(path: &Path) -> Result<Vec<Sam>, String> {
	let file = match File::open(&path) {
		Err(why) => return Err(format!("open {}: {}", path.display(), why)),
		Ok(file) => file,
	};
	let reader = BufReader::new(&file);
	reader.lines().try_fold(Vec::new(), |mut sam, line| {
		let line = line.unwrap();
		let field: Vec<_> = line.split('\t').collect();

		// Ignore the header section
		if field[0].chars().next().unwrap() != '@' {
			let flag: u16 = field[1].parse().unwrap();
			// Primary alignment
			if flag == 0 || flag == 16 {
				sam.push(Sam::Mapped(Mapped {
					qname: field[0].to_string(),
					len: field[9].len() as u32,
					strand: (flag == 16),
					rname: field[2].to_string(),
					pos_min: field[3].parse().unwrap(),
					pos_max: field[3].parse().unwrap(),
					mapq: field[4].parse().unwrap(),
					als: {
						let al_s: Vec<_> = field[13].split(':').collect();
						al_s[2].parse().unwrap()
					},
					secondaries: Vec::new(),
				}));
			}
			// Secondary alignment
			else if (flag & 256) != 0 {
				let len = sam.len();
				let entry = sam.get_mut(len - 1).unwrap();
				match entry {
					Sam::Mapped(x) => x.secondaries.push(Secondary {
						strand: (flag & 16) == 1,
						rname: field[2].to_string(),
						pos: field[3].parse().unwrap(),
						als: {
							let al_s: Vec<_> = field[13].split(':').collect();
							al_s[2].parse().unwrap()
						},
					}),
					Sam::Unmapped(x) => {
						return Err(format!(
							"Unmapped sequence with a secondary alignment: {} {}",
							path.display(),
							x.qname
						));
					}
				}
			}
			// Supplementary alignment
			else if (flag & 2048) != 0 {
				let pos = field[3].parse().unwrap();
				let len = sam.len();
				let entry = sam.get_mut(len - 1).unwrap();
				match entry {
					Sam::Mapped(x) => {
						if x.pos_min > pos {
							x.pos_min = pos;
						} else if x.pos_max < pos {
							x.pos_max = pos;
						}
					}
					Sam::Unmapped(x) => {
						return Err(format!(
							"Unmapped sequence with a supplementary alignment: {} {}",
							path.display(),
							x.qname
						));
					}
				}
			}
			// Unmapped
			else if flag == 4 {
				sam.push(Sam::Unmapped(Unmapped {
					qname: field[0].to_string(),
					len: field[9].len() as u32,
				}))
			} else {
				return Err(format!("Unknown flag: {} {}", path.display(), line));
			}
		}
		Ok(sam)
	})
}

fn compare_sam(
	tgt: &Vec<Sam>,
	test: &Vec<Sam>,
	distance: f32,
	qualities: &Vec<u8>,
	output: &Option<String>,
	p_only: bool,
) -> Result<(), String> {
	let mut iter = zip(tgt, test);

	let mut files: Option<(File, File, File)> = match output {
		Some(output) => {
			let name = format!("{}_gain.txt", output);
			let path = Path::new(&name);
			let gain_file = match File::create(path) {
				Err(err) => return Err(format!("{}", err)),
				Ok(f) => f,
			};

			let name = format!("{}_loss.txt", output);
			let path = Path::new(&name);

			let loss_file = match File::create(path) {
				Err(err) => return Err(format!("{}", err)),
				Ok(f) => f,
			};

			let name = format!("{}_diff.txt", output);
			let path = Path::new(&name);
			let diff_file = match File::create(path) {
				Err(err) => return Err(format!("{}", err)),
				Ok(f) => f,
			};

			Some((gain_file, loss_file, diff_file))
		}
		None => None,
	};

	let qualities_len = qualities.len();

	let mut gain: Vec<u32> = vec![0; qualities_len];
	let mut loss: Vec<u32> = vec![0; qualities_len];
	let mut diff: Vec<u32> = vec![0; qualities_len];

	iter.try_for_each(|x| {
		let tgt = x.0;
		let test = x.1;
		match tgt {
			Sam::Mapped(tgt) => match test {
				Sam::Mapped(test) => {
					let distance = (distance * tgt.len as f32).ceil() as u32;
					if test.rname != tgt.rname
						|| test.strand != tgt.strand
						|| test.pos_max < tgt.pos_min - distance
						|| test.pos_min > tgt.pos_max + distance
					{
						let res = test.secondaries.iter().try_for_each(|x| {
							if x.rname == tgt.rname {
								None
							} else {
								Some(())
							}
						});
						match res {
							None => {}
							Some(_) => {}
						}
					}
				}
				Sam::Unmapped(_) => {
					increase_counter(&mut loss, &qualities, tgt.mapq);
					if let Some((_, file, _)) = &mut files {
						match writeln!(file, "{}", tgt) {
							Err(err) => return Err(format!("{}", err)),
							Ok(_) => {}
						};
					}
				}
			},
			Sam::Unmapped(_) => match test {
				Sam::Mapped(test) => {
					increase_counter(&mut gain, &qualities, test.mapq);
					if let Some((file, _, _)) = &mut files {
						match writeln!(file, "{}", tgt) {
							Err(err) => return Err(format!("{}", err)),
							Ok(_) => {}
						};
					}
				}
				Sam::Unmapped(_) => {}
			},
		}
		Ok(())
	})?;

	Ok(())
}

#[inline]
fn increase_counter(count: &mut Vec<u32>, qualities: &Vec<u8>, quality: u8) {
	let iter = zip(count, qualities);
	for i in iter {
		if quality > *i.1 {
			*i.0 += 1;
			break;
		}
	}
}
