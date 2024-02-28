#!/bin/bash

make clean
mkdir -p kat

for iut in \
	PLOVER_128_1	PLOVER_128_2	PLOVER_128_4	\
	PLOVER_128_8	PLOVER_128_16	PLOVER_128_32	\
	PLOVER_192_1	PLOVER_192_2	PLOVER_192_4	\
	PLOVER_192_8	PLOVER_192_16	PLOVER_192_32	\
	PLOVER_256_1	PLOVER_256_2	PLOVER_256_4	\
	PLOVER_256_8	PLOVER_256_16	PLOVER_256_32
do
	echo	=== $iut ===
	make -f Makefile.kat IUT="$iut" obj-clean
	make -f Makefile.kat IUT="$iut"
	cd kat
	./xgen_$iut
	cd ..
done
