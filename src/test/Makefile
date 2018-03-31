.PHONY: default all

default: all

falcon: ../c/DW_banded.c ../c/kmer_lookup.c ../c/falcon.c falcon_main.c  ../c/common.h
	gcc  ../c/DW_banded.c ../c/kmer_lookup.c ../c/falcon.c  falcon_main.c -I../c/ -o falcon

test.out: falcon
	cat test.in | ./falcon > test.out
	diff test.out test.expected && echo test passed

test.in: generate_test.py
	python3 generate_test.py > test.in

all: test.out

clean:
	rm test.out falcon
