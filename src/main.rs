extern crate getopts;
use getopts::Options;
use std::env::args;
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::iter::zip;
use std::path::Path;
use std::process;
use std::str::FromStr;

//TODO: check if there are also supplementatry alignments for secondaries
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
		"Comparison mode [all]\n\
		- all: match the primary and the secondary alignments of the test file with the primary and the secondary of the target file respectively\n\
		- prim_tgt: match the primary and the secondary alignments of the test file with the primar aligments of the target file\n\
		- prim: match the primary alignments of the test file with the primary aligments of the target file",
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

	let distance: f32 = match matches.opt_get_default("d", 1.0) {
		Ok(d) => d,
		Err(err) => {
			eprintln!("[ERROR]: {} 'd'", err);
			process::exit(1);
		}
	};

	let qualities: Vec<u8> = match matches.opt_str("q") {
		Some(q) => q
			.split(',')
			.map(|x| match x.parse() {
				Ok(u) => u,
				Err(err) => {
					eprintln!("[ERROR]: {} 'q'", err);
					process::exit(1);
				}
			})
			.collect(),
		None => [60, 10, 1, 0].to_vec(),
	};

	let mode = match matches.opt_get_default("m", Mode::All) {
		Ok(m) => m,
		Err(err) => {
			eprintln!("[ERROR]: {}", err);
			process::exit(1);
		}
	};

	if matches.free.len() != 2 {
		print_usage(&program, opts);
		return;
	}

	let path = Path::new(&matches.free[0]);
	let sam_tgt = match parse_sam(&path) {
		Ok(sam_tgt) => sam_tgt,
		Err(err) => {
			eprintln!("[ERROR]: {} {}", path.display(), err);
			process::exit(1);
		}
	};

	let path = Path::new(&matches.free[1]);
	let sam_test = match parse_sam(&path) {
		Ok(sam_test) => sam_test,
		Err(err) => {
			eprintln!("[ERROR]: {} {}", path.display(), err);
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

	compare_sam(&sam_tgt, &sam_test, distance, &qualities, &output, mode);
}

fn print_usage(program: &str, opts: Options) {
	let brief = format!("Usage: {} [options] <target.sam> <test.sam>", program);
	print!("{}", opts.usage(&brief));
}

#[derive(Copy, Clone)]
enum Mode {
	All,
	PrimTgt,
	Prim,
}

impl FromStr for Mode {
	type Err = String;

	fn from_str(s: &str) -> Result<Self, Self::Err> {
		if s == "all" {
			Ok(Mode::All)
		} else if s == "prim_tgt" {
			Ok(Mode::PrimTgt)
		} else if s == "prim" {
			Ok(Mode::Prim)
		} else {
			Err("Unrecognized mode 'm'".to_string())
		}
	}
}

#[derive(Clone)]
struct Secondary {
	rname: String,
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
			"{}\t{}\t{}\t{}\t{}\t{}\t{}",
			self.qname,
			self.len,
			self.rname,
			self.pos_min,
			self.pos_max,
			if self.strand { '-' } else { '+' },
			self.mapq,
		)?;
		self.secondaries.iter().try_for_each(|x| {
			write!(
				f,
				"\t[{}\t{}\t{}\t{}]",
				x.rname,
				x.pos,
				if x.strand { '-' } else { '+' },
				x.als
			)
		})?;
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
		Err(why) => return Err(format!("open: {}", why)),
		Ok(file) => file,
	};
	let reader = BufReader::new(&file);
	reader.lines().try_fold(Vec::new(), |mut sam, line| {
		let line = match line {
			Err(err) => return Err(format!("line: {}", err)),
			Ok(l) => l,
		};
		let field: Vec<_> = line.split('\t').collect();

		// Ignore the header section
		if field[0].chars().next().ok_or("parse_sam")? != '@' {
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
	mode: Mode,
) {
	let iter = zip(tgt, test);

	let mut files: Option<(File, File, File)> = match output {
		Some(output) => {
			let name = format!("{}_gain.txt", output);
			let path = Path::new(&name);
			let gain_file = match File::create(path) {
				Err(err) => {
					eprintln!("[ERROR]: create {}", err);
					process::exit(1);
				}
				Ok(f) => f,
			};

			let name = format!("{}_loss.txt", output);
			let path = Path::new(&name);

			let loss_file = match File::create(path) {
				Err(err) => {
					eprintln!("[ERROR]: create {}", err);
					process::exit(1);
				}
				Ok(f) => f,
			};

			let name = format!("{}_diff.txt", output);
			let path = Path::new(&name);
			let diff_file = match File::create(path) {
				Err(err) => {
					eprintln!("[ERROR]: create {}", err);
					process::exit(1);
				}
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

	iter.for_each(|x| {
		let tgt = x.0;
		let test = x.1;
		match tgt {
			Sam::Mapped(tgt) => match test {
				Sam::Mapped(test) => {
					let distance = (distance * tgt.len as f32).ceil() as u32;
					let file = match &mut files {
						Some((_, _, file)) => Some(file),
						None => None,
					};

					match mode {
						Mode::All => {
							compare_all(&tgt, &test, &mut diff, &qualities, file, distance)
						}
						Mode::PrimTgt => {
							compare_prim_tgt(&tgt, &test, &mut diff, &qualities, file, distance)
						}
						Mode::Prim => {
							compare_prim(&tgt, &test, &mut diff, &qualities, file, distance)
						}
					}
				}
				Sam::Unmapped(_) => {
					increase_counter(&mut loss, &qualities, tgt.mapq);
					if let Some((_, file, _)) = &mut files {
						if let Err(err) = writeln!(file, "{}", tgt) {
							eprintln!("[ERROR]: write {}", err);
							process::exit(1);
						}
					}
				}
			},
			Sam::Unmapped(_) => match test {
				Sam::Mapped(test) => {
					increase_counter(&mut gain, &qualities, test.mapq);
					if let Some((file, _, _)) = &mut files {
						if let Err(err) = writeln!(file, "{}", tgt) {
							eprintln!("[ERROR]: write {}", err);
							process::exit(1);
						}
					}
				}
				Sam::Unmapped(_) => {}
			},
		}
	});

	let iter = zip(qualities, gain);
	iter.for_each(|x| println!("G\t{}\t{}", x.0, x.1));
	println!("");
	let iter = zip(qualities, loss);
	iter.for_each(|x| println!("L\t{}\t{}", x.0, x.1));
	let iter = zip(qualities, diff);
	println!("");
	iter.for_each(|x| println!("D\t{}\t{}", x.0, x.1));
}

#[inline]
fn increase_counter(count: &mut Vec<u32>, qualities: &Vec<u8>, quality: u8) {
	let iter = zip(count, qualities);
	for i in iter {
		if quality >= *i.1 {
			*i.0 += 1;
			break;
		}
	}
}
fn compare_all(
	tgt: &Mapped,
	test: &Mapped,
	count: &mut Vec<u32>,
	qualities: &Vec<u8>,
	file: Option<&mut File>,
	distance: u32,
) {
	if test.rname != tgt.rname
		|| test.strand != tgt.strand
		|| test.pos_max + distance < tgt.pos_min
		|| test.pos_min > tgt.pos_max + distance
	{
		increase_counter(count, qualities, tgt.mapq);
		if let Some(file) = file {
			if let Err(err) = writeln!(file, ">>>>>>>>\n{}\n<<<<<<<<\n{}", tgt, test) {
				eprintln!("[ERROR]: write {}", err);
				process::exit(1);
			}
		}
		return;
	}

	for tgt_s in &tgt.secondaries {
		let mut absent = true;
		for test_s in &test.secondaries {
			if tgt_s.rname == test_s.rname
				&& tgt_s.strand == test_s.strand
				&& test_s.pos + distance >= tgt_s.pos
				&& test_s.pos <= tgt_s.pos + distance
			{
				absent = false;
				break;
			}
		}
		if absent {
			increase_counter(count, qualities, tgt.mapq);
			if let Some(file) = file {
				if let Err(err) = writeln!(file, ">>>>>>>>\n{}\n<<<<<<<<\n{}", tgt, test) {
					eprintln!("[ERROR]: write {}", err);
					process::exit(1);
				}
			}
			return;
		}
	}
}

fn compare_prim_tgt(
	tgt: &Mapped,
	test: &Mapped,
	count: &mut Vec<u32>,
	qualities: &Vec<u8>,
	file: Option<&mut File>,
	distance: u32,
) {
	if test.rname != tgt.rname
		|| test.strand != tgt.strand
		|| test.pos_max + distance < tgt.pos_min
		|| test.pos_min > tgt.pos_max + distance
	{
		test.secondaries.iter().for_each(|s| {
			if s.rname == tgt.rname
				&& s.strand == tgt.strand
				&& s.pos >= tgt.pos_min - distance
				&& s.pos <= tgt.pos_max + distance
			{
				return;
			}
		});
		increase_counter(count, qualities, tgt.mapq);
		if let Some(file) = file {
			if let Err(err) = writeln!(file, ">>>>>>>>\n{}\n<<<<<<<<\n{}", tgt, test) {
				eprintln!("[ERROR]: write {}", err);
				process::exit(1);
			}
		}
	}
}

fn compare_prim(
	tgt: &Mapped,
	test: &Mapped,
	count: &mut Vec<u32>,
	qualities: &Vec<u8>,
	file: Option<&mut File>,
	distance: u32,
) {
	if test.rname != tgt.rname
		|| test.strand != tgt.strand
		|| test.pos_max + distance < tgt.pos_min
		|| test.pos_min > tgt.pos_max + distance
	{
		increase_counter(count, qualities, tgt.mapq);
		if let Some(file) = file {
			if let Err(err) = writeln!(file, ">>>>>>>>\n{}\n<<<<<<<<\n{}", tgt, test) {
				eprintln!("[ERROR]: write {}", err);
				process::exit(1);
			}
		}
	}
}
