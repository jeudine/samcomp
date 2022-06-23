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
		"Generate gain, loss and diff files and output them with the prefix NAME",
		"NAME",
	);
	let matches = match opts.parse(&args[1..]) {
		Ok(m) => m,
		Err(f) => panic!("{}", f.to_string()),
	};

	if matches.opt_present("h") {
		print_usage(&program, opts);
		return;
	}

	let output = matches.opt_str("o");

	if matches.free.len() != 2 {
		print_usage(&program, opts);
		return;
	}

	let path = Path::new(&matches.free[0]);
	let sam_tgt = parse_sam(&path);
	let path = Path::new(&matches.free[1]);
	let sam_test = parse_sam(&path);

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

	match output {
		Some(output) => compare_sam_output(&sam_tgt, &sam_test, &output),
		None => compare_sam(&sam_tgt, &sam_test),
	}

	/*
	let (same, gain, loss, diff) = iter.fold((0, 0, 0, 0), |(same, gain, loss, diff), ref x| {
		if x.0.flag == 4 {
			if x.1.flag == 4 {
				(same, gain, loss, diff)
			} else {
				writeln!(&mut loss_file, "{}", x.1.qname).unwrap();
				(same, gain, loss + 1, diff)
			}
		} else {
			if x.1.flag == 4 {
				(same, gain + 1, loss, diff)
			} else {
				if x.0.rname.eq(&x.1.rname)
					&& x.0.flag == x.1.flag
					&& (x.1.pos as i32 - x.0.pos as i32).abs() <= 10000
				{
					(same + 1, gain, loss, diff)
				} else {
					writeln!(&mut diff_file, "{}", x.1.qname).unwrap();
					(same, gain, loss, diff + 1)
				}
			}
		}
	});
	let total = same + gain + loss + diff;
	println!("same: {}%", (same * 100) as f32 / total as f32);
	println!("gain: {}%", (gain * 100) as f32 / total as f32);
	println!("loss: {}%", (loss * 100) as f32 / total as f32);
	println!("different: {}%", (diff * 100) as f32 / total as f32);
	*/
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
		);
		self.secondaries.iter().for_each(|x| {
			write!(f, "\t[{}\t{}\t{}\t{}]", x.rname, x.pos, x.strand, x.als);
		});
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

fn parse_sam(path: &Path) -> Vec<Sam> {
	let file = match File::open(&path) {
		Err(why) => panic!("open {}: {}", path.display(), why),
		Ok(file) => file,
	};
	let reader = BufReader::new(&file);
	reader.lines().fold(Vec::new(), |mut sam, line| {
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
						panic!("Unmapped sequence with a secondary alignment: {}", x.qname)
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
						panic!(
							"Unmapped sequence with a supplementary alignment: {}",
							x.qname
						)
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
				panic!("Unknown flag: {}", line);
			}
		}
		sam
	})
}

fn compare_sam_output(tgt: &Vec<Sam>, test: &Vec<Sam>, output: &str) {
	let name = format!("{}_gain.txt", output);
	let path = Path::new(&name);
	let mut gain_file = File::create(path).unwrap();

	let name = format!("{}_loss.txt", output);
	let path = Path::new(&name);
	let mut loss_file = File::create(path).unwrap();

	let name = format!("{}_diff.txt", output);
	let path = Path::new(&name);
	let mut diff_file = File::create(path).unwrap();

	let iter = zip(tgt, test);

	iter.for_each(|x| {
		let tgt = x.0;
		let test = x.1;

		match tgt {
			Sam::Mapped(tgt) => match test {
				Sam::Mapped(test) => {
					if test.rname != tgt.rname
						|| test.strand != tgt.strand
						|| test.pos_max < tgt.pos_min - tgt.len
						|| test.pos_min > tgt.pos_max + tgt.len
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
				Sam::Unmapped(_) => writeln!(&mut loss_file, "{}", tgt).unwrap(),
			},
			Sam::Unmapped(_) => match test {
				Sam::Mapped(test) => writeln!(&mut gain_file, "{}", test).unwrap(),
				Sam::Unmapped(_) => {}
			},
		}
	});
}

fn compare_sam(tgt: &Vec<Sam>, test: &Vec<Sam>) {
	let iter = zip(tgt, test);

	iter.for_each(|x| {
		let tgt = x.0;
		let test = x.1;

		match tgt {
			Sam::Mapped(tgt) => match test {
				Sam::Mapped(test) => {
					if test.rname != tgt.rname
						|| test.strand != tgt.strand
						|| test.pos_max < tgt.pos_min - tgt.len
						|| test.pos_min > tgt.pos_max + tgt.len
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
				Sam::Unmapped(_) => {}
			},
			Sam::Unmapped(_) => match test {
				Sam::Mapped(test) => {}
				Sam::Unmapped(_) => {}
			},
		}
	});
}
