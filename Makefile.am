bin_PROGRAMS=pangenomesim
pangenomesim_SOURCES=src/global.h src/pangenomesim.cxx src/util.h src/util.cxx src/img.h src/img.cxx src/locus.h src/locus.cxx
pangenomesim_CPPFLAGS=-Wall -Wextra -Weffc++ -std=c++14
pangenomesim_CXXFLAGS=

clean-local:
	rm -rf *.fasta *.maf *.seed
	rm -rf foo

.PHONY: format
format:
	clang-format -i src/*.cxx src/*.h
