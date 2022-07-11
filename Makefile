all:
	cargo build --release

test:
	#cargo run --release -- test/target.sam test/test.sam -d 0.5 -q60,30,20,10,5,2,1,0
	cargo run --release -- test/target.sam test/test.sam -q60,30,20,10,5,2,1,0 -m all -o test -x

clean:
	cargo clean

.PHONY: all test clean
