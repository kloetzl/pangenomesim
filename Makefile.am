bin_PROGRAMS=pangenomesim
pangenomesim_SOURCES=src/global.h src/pangene.cxx src/pangene.h src/pangenomesim.cxx src/simple.h src/simple.cxx src/util.h src/util.cxx src/evo_model.h src/evo_model.cxx src/img.h src/img.cxx src/locus.h src/locus.cxx
pangenomesim_CPPFLAGS=-std=c++14 -Wall -Wextra -Weffc++
pangenomesim_CXXFLAGS=-O1

clean-local:
	rm -rf *.fasta *.maf *.seed
	rm -rf foo

.PHONY: format
format:
	clang-format -i src/*.cxx src/*.h
