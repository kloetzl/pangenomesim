#include "global.h"
#include "img.h"
#include "locus.h"
#include "util.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <regex>
#include <string>
#include <vector>

extern "C" {
#include <err.h>
#include <errno.h>
#include <getopt.h>
}

#ifndef VERSION
#define VERSION "v0.2"
#endif

void usage(int);
void version();

std::string OUT_DIR = std::string("./");
std::default_random_engine RNG;
bool VERBOSE = false;

int main(int argc, char *argv[])
{
	static const struct option long_options[] = {
		{"version", no_argument, NULL, 0},
		{"help", no_argument, NULL, 0},
		{"param", required_argument, NULL, 'p'},
		{"out-dir", required_argument, NULL, 'o'},
		{"verbose", no_argument, NULL, 'v'},
		{0, 0, 0, 0} // no comment
	};

	auto model = img_model();

	while (true) {
		int option_index = 0;
		int c = getopt_long(argc, argv, "o:p:v", long_options, &option_index);

		if (c == -1) {
			break;
		}

		switch (c) {
			case 0: {
				auto option_string =
					std::string(long_options[option_index].name);
				if (option_string == "help") {
					usage(EXIT_SUCCESS);
				} else if (option_string == "version") {
					version();
				} else {
					return 1;
				}
				break;
			}
			case 'o': {
				OUT_DIR = std::string(optarg) + "/";
				break;
			}
			case 'p': {
				auto arg_string = std::string(optarg);

				// implement `getsubsopt(3)` like functionality
				auto pattern = std::string("^(\\w+)(?:=([a-z_0-9.]+))?,?");
				auto r = std::regex(pattern);

				auto m = std::smatch();
				while (std::regex_search(arg_string, m, r)) {
					auto key = m[1];
					auto value = m.size() == 3 ? m[2] : std::string();

					// may throw
					model.parse_param(key, value);
					arg_string = m.suffix();
				}

				if (!arg_string.empty()) {
					errx(1, "invalid parameter %s", optarg);
				}

				break;
			}
			case 'v': {
				VERBOSE = true;
				break;
			}
			case '?': /* intentional fall-thorugh */
			default: usage(EXIT_FAILURE); break;
		}
	}

	// create ouput directory
	mkpath(OUT_DIR);

	model.simulate(); // do the work

	{ // reproducible information
		auto rep_file_name = OUT_DIR + "reproducible.seed";
		auto rep_file = std::ofstream(rep_file_name);
		check_io(rep_file, rep_file_name);

		rep_file << "pangenomesim " VERSION << "\n\n";
		rep_file << "full options:\n";
		rep_file << model.parameters() << std::endl;

		check_io(rep_file, rep_file_name);
	}

	{ // print reference, core-, and accessory-genome
		auto ref_file_name = OUT_DIR + "ref.fasta";
		auto cor_file_name = OUT_DIR + "core.fasta";
		auto acc_file_name = OUT_DIR + "accessory.fasta";
		auto ref_fasta_file = std::ofstream(ref_file_name);
		auto cor_fasta_file = std::ofstream(cor_file_name);
		auto acc_fasta_file = std::ofstream(acc_file_name);
		check_io(ref_fasta_file, ref_file_name);
		check_io(cor_fasta_file, cor_file_name);
		check_io(acc_fasta_file, acc_file_name);

		for (const auto &loc : model.get_reference()) {
			ref_fasta_file << loc.to_fasta();
		}

		for (const auto &loc : model.get_core()) {
			cor_fasta_file << loc.to_fasta();
		}

		for (const auto &loc : model.get_accessory()) {
			acc_fasta_file << loc.to_fasta();
		}

		check_io(ref_fasta_file, ref_file_name);
		check_io(cor_fasta_file, cor_file_name);
		check_io(acc_fasta_file, acc_file_name);
	}

	// print derived genomes
	for (size_t gn = 0; gn < model.get_num_genomes(); gn++) {
		auto file_name = OUT_DIR + genome_name(gn) + ".fasta";
		auto fasta_file = std::ofstream(file_name);
		check_io(fasta_file, file_name);

		auto loci = model.get_genome(gn);
		for (const auto &loc : loci) {
			fasta_file << loc.to_fasta();
		}

		check_io(fasta_file, file_name);
	}

	{ // compute gene frequency spectrum
		auto num_loci = model.get_num_loci();
		auto gfs = std::vector<size_t>(model.get_num_genomes() + 1);
		for (size_t i = 0; i < num_loci; i++) {
			auto that = model.get_locus(i);
			gfs.at(that.size())++;
		}

		auto file_name = OUT_DIR + "genefrequency.gfs";
		auto gfs_file = std::ofstream(file_name);
		check_io(gfs_file, file_name);

		for (auto it = gfs.begin() + 1; it != gfs.end(); it++) {
			gfs_file << *it << " ";
		}
		gfs_file << std::endl;

		check_io(gfs_file, file_name);
	}

	{ // coalescent
		auto file_name = OUT_DIR + "coalescent.newick";
		auto tree_file = std::ofstream(file_name);
		check_io(tree_file, file_name);

		auto pre = [&tree_file](const tree_node &self) {
			if (self.is_branch()) {
				tree_file << "(";
			}
		};
		auto process = [&tree_file](const tree_node &self) {
			if (self.is_leaf()) {
				tree_file << genome_name(self.get_index());
			} else {
				tree_file << ",";
			}
		};
		auto post = [&tree_file](const tree_node &self) {
			if (self.is_branch()) {
				tree_file << ")";
			}
			if (self.get_time() != 0.0) {
				tree_file << ":" << self.get_time();
			}
		};

		model.get_root().traverse(pre, process, post);

		tree_file << ";\n";

		check_io(tree_file, file_name);
	}

	{ // matrix
		auto file_name = OUT_DIR + "distance.mat";
		auto mat_file = std::ofstream(file_name);
		check_io(mat_file, file_name);

		auto n = model.get_num_genomes();
		const auto &distmatrix = model.get_distmatrix();

		assert(distmatrix.size() == n);

		mat_file << n << "\n";
		for (size_t i = 0; i < n; i++) {
			mat_file << genome_name(i);
			for (auto d : distmatrix[i]) {
				mat_file << "  " << std::setw(8) << std::setprecision(4) << d;
			}
			mat_file << "\n";
		}

		check_io(mat_file, file_name);
	}

	{ // MAF alignment
		auto file_name = OUT_DIR + "alignment.maf";
		auto maf_file = std::ofstream(file_name);
		check_io(maf_file, file_name);
		maf_file << "##maf version=1 program=pangenomesim\n";

		auto loci_length = model.get_loci_length();

		auto num_genomes = model.get_num_genomes();
		// compute sequence sizes
		auto seq_sizes = std::vector<size_t>(num_genomes);
		for (size_t i = 0; i < num_genomes; i++) {
			seq_sizes[i] = model.get_genome(i).size() * loci_length;
		}

		auto ref = model.get_reference();
		auto ref_size = ref.size() * loci_length;

		auto s_line = [&loci_length, &maf_file, &seq_sizes](const auto &loc) {
			auto gid = loc.get_genome_id();
			maf_file << "s " << genome_name(gid)  // genome ID
					 << "." << loc.get_locus_id() // contig ID
					 << " 0"					  // start
					 << " " << loci_length		  // size
					 << " +"					  // strand
					 << " "
					 << seq_sizes[gid] // size of the entire source sequence
					 << " " << loc.get_nucl() // sequence
					 << "\n";
		};

		for (size_t i = 0; i < model.get_num_loci(); i++) {
			maf_file << "a\n";
			// print reference
			maf_file << "s " << genome_name(-1)				 // genome ID
					 << "." << ref[i].get_locus_id() << " 0" // start
					 << " " << loci_length					 // size
					 << " +"								 // strand
					 << " " << ref_size						 // size s.a.
					 << " " << ref[i].get_nucl()			 // sequence
					 << "\n";

			auto tmp = model.get_locus(i);
			std::for_each(tmp.begin(), tmp.end(), s_line);

			maf_file << std::endl;
		}

		check_io(maf_file, file_name);
	}

	return 0;
}

void usage(int exit_code)
{
	static const char str[] = {
		"usage: pangenomesim [OPTIONS...]\n"
		"Simulate a pan-genome.\n\n"
		"  -o, --out-dir DIR        The directory to write files to\n"
		"  -p, --param KEY=VALUE    Set simulation parameter\n"
		"  -v, --verbose            Print additional information\n"
		"      --help               Display help and exit\n"
		"      --version            Output version information and exit\n" //
	};

	fprintf(exit_code == EXIT_SUCCESS ? stdout : stderr, "%s", str);
	exit(exit_code);
}

void version()
{
	static const char str[] = {
		"pangenomesim " VERSION "\n"
		"Copyright (C) 2017,  Fabian Klötzl "
		"<fabian-pangenomesim@kloetzl.info>\n"
		"BSD 2-Clause License\n" //
	};
	printf("%s", str);
	exit(EXIT_SUCCESS);
}
