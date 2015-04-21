DW_align.so: DW_banded.c common.h
	gcc DW_banded.c -O3 -shared -fPIC -o DW_align.so

kmer_lookup.so: kmer_lookup.c common.h
	gcc kmer_lookup.c -O3 -shared -fPIC -o kmer_lookup.so

#falcon: DW_banded.c common.h kmer_lookup.c falcon.c 
#	gcc DW_banded.c kmer_lookup.c falcon.c -O4 -o falcon -fPIC 

falcon.so: falcon.c common.h DW_banded.c kmer_lookup.c
	gcc DW_banded.c kmer_lookup.c falcon.c -O3 -shared -fPIC -o falcon.so 

#falcon2.so: falcon.c common.h DW_banded_2.c kmer_lookup.c
#	gcc DW_banded_2.c kmer_lookup.c falcon.c -O3 -shared -fPIC -o falcon2.so 

clean:
	rm  *.so

all: DW_align.so kmer_lookup.so falcon.so

