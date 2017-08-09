bin_PROGRAMS=pangenomesim
pangenomesim_SOURCES=src/global.h src/pangenomesim.cxx src/util.h src/util.cxx src/img.h src/img.cxx src/gene.h src/gene.cxx
pangenomesim_CPPFLAGS=-Wall -Wextra -Weffc++ -std=c++14
pangenomesim_CXXFLAGS=

dist_man_MANS=pangenomesim.1

clean-local:
	rm -rf *.fasta *.maf *.seed
	rm -rf foo

.PHONY: format
format:
	clang-format -i src/*.cxx src/*.h
