all:
	cargo build --release

test:
	cargo run --release -- test/target.sam test/target.sam

clean:
	cargo clean

.PHONY: all test clean
