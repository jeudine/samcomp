# samcomp

A comparison tool for Sequence Alignment/Map files written in Rust.

## Installation

### From source
if you want to build samcomp from source, you need Rust. You can then use `cargo` to build everything:

```bash
cargo install samcomp
```

## Usage

```bash
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
                        - all: match the primary and the secondary alignments
                        of the test file with the primary and the secondary of
                        the target file respectively
                        - prim_tgt: match the primary and the secondary
                        alignments of the test file with the primar aligments
                        of the target file
                        - prim: match the primary alignments of the test file
                        with the primary aligments of the target file
```
