bin_PROGRAMS=pangenomesim
pangenomesim_SOURCES=src/global.h src/pangenomesim.cxx src/util.h src/util.cxx src/img.h src/img.cxx src/gene.h src/gene.cxx src/pfasta.h src/pfasta.c
pangenomesim_CPPFLAGS=-Wall -Wextra
pangenomesim_CXXFLAGS= -std=c++14
pangenomesim_CFLAGS=

dist_man_MANS=pangenomesim.1

clean-local:
	rm -rf *.fasta *.maf *.seed
	rm -rf foo

.PHONY: format
format:
	clang-format -i src/*.cxx src/*.h
