# samcomp

A comparison tool for Sequence Alignment/Map files written in Rust.

## Installation

### From source
if you want to build samcomp from source, you need Rust. You can then use `cargo` to build everything:

```bash
cargo install samcomp
```

## Usage

```
Usage: samcomp [options] <target.sam> <test.sam>

Options:
    -h, --help          print this help menu
    -o NAME             generate gain, loss and diff files with the name of
                        the reads and output them with the prefix NAME
    -d FLOAT            a location in the tested file is considered to be
                        similar to the one in the target file if the distance
                        between both locations is less than FLOAT fraction of
                        the read length [1.0]
    -q UINT1,UINT2,...  output the results using UINT1,UINT2,... (such as
                        UINT1 > UINT2 > ...) as the quality thresholds
                        [60,10,1,0]
    -m STR              Comparison mode [all]
                        - all: compare the primary and the secondary
                        alignments of the tested file with the primary and the
                        secondary alignments of the target file respectively
                        - prim_tgt: compare the primary, the secondary and the
                        supplementary alignments of the tested file with the
                        primary aligments of the target file
                        - prim: compare the primary alignments of the tested
                        file with the primary aligments of the target file
                        - prim_supp: compare the primary and the supplementary
                        alignments of the tested file with the primary
                        aligments of the target file
```

samcomp evaluates the differences between 2 SAM files (target file and tested file) containing the same reads.

Here is an example output:

```
M	60	6	3
M	30	7	3
M	20	7	3
M	10	7	3
M	5	7	3
M	2	7	9
M	1	8	9
M	0	8	9

G	60	0
G	30	0
G	20	0
G	10	0
G	5	0
G	2	1
G	1	1
G	0	1

L	60	0
L	30	0
L	20	0
L	10	0
L	5	0
L	2	0
L	1	0
L	0	0

D	60	1
D	30	1
D	20	1
D	10	1
D	5	1
D	2	1
D	1	2
D	0	2
```

Each M-line (Mapped) gives the number of mapped reads in the target file (col 2) and the tested file (col 3) with a mapping quality equal to or greater than the threshold (col 1).

Each G-line (Gain) gives the number of reads (col 2) which are mapped in the tested file and unmapped in the target file with a mapping quality equal to or greater than the threshold (col 1).

Each L-line (Loss) gives the number of reads (col 2) which are unmapped in the tested file and mapped in the target file with a mapping quality equal to or greater than the threshold (col 1).

Each D-line (Difference) gives the number of reads (col 2) which have a different mapping location in the target and tested file with a target mapping quality equal to or greater than the threshold (col 1).

***

With the option `-o`, it is also possible to output the names (one per line) of the reads belonging to the G, L, and D categories.
