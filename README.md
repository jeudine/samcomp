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
    -o NAME             generate gain, loss and diff files and output them
                        with the prefix NAME
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
                        - prim_tgt: compare the primary and the secondary
                        alignments of the tested file with the primary
                        aligments of the target file
                        - prim: compare the primary alignments of the tested
                        file with the primary aligments of the target file
```

samcomp evaluates the differences between 2 SAM files (target file and tested file) containing the same reads.

Here is an example output:

```
G	60	5
G	10	0
G	1	69
G	0	0

L	60	0
L	10	7
L	1	0
L	0	0

D	60	106
D	10	70
D	1	258
D	0	154
```

Each G-line (Gain) gives the number of reads (row 2) which are mapped in the tested file and unmapped in the target file with a mapping quality equal to or greater than the threshold (row 1).

Each L-line (Loss) gives the number of reads (row 2) which are unmapped in the tested file and mapped in the target file with a mapping quality equal to or greater than the threshold (row 1).

Each D-line (Difference) gives the number of reads (row 2) which have a different mapping location in the target and tested file with a target mapping quality equal to or greater than the threshold (row 1).
