bin_PROGRAMS=pangenomesim
pangenomesim_SOURCES=src/global.h src/pangene.cxx src/pangene.h src/pangenomesim.cxx src/simulate.h src/simulate.cxx
pangenomesim_CPPFLAGS=-std=c++14 -Wall -Wextra -Weffc++
pangenomesim_CXXFLAGS=

clean-local:
	rm -rf *.fasta *.maf *.seed
	rm -rf foo

.PHONY: format
format:
	clang-format -i src/*.cxx src/*.h
