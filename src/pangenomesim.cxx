#include "global.h"
#include "pangene.h"
#include <algorithm>
#include <cstdio>
#include <err.h>
#include <errno.h>
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <random>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

#ifndef VERSION
#define VERSION "v0.1"
#endif

void usage(int);
void version();

double EVO_DISTANCE = 0.01;
size_t GENE_LENGTH = 100;
size_t NUM_GENOMES = 3;
size_t NUM_GENES = 5;
std::string OUT_DIR = std::string("./");
double PROBABILITY = 0.8;
size_t SEED = 1729;
double THRESHOLD = 1.0;
std::default_random_engine RNG;
bool VERBOSE = false;

auto create_distribution()
{
	auto distribution = std::vector<double>(NUM_GENOMES + 1);

	auto tail = std::pow(1 - PROBABILITY, NUM_GENOMES);
	auto factor = 1 / (1 - tail);
	std::cout << "factor\t" << factor << "\n";

	for (size_t i = 1; i < NUM_GENOMES + 1; i++) {
		distribution[i] =
			distribution[i - 1] +
			PROBABILITY * std::pow(1 - PROBABILITY, i - 1) * factor;
	}

	return distribution;
}

auto genome_name(size_t genome_number)
{
	auto ret = std::string("genome");
	auto num = std::to_string(genome_number);
	auto max_digits = std::ceil(std::log10(NUM_GENOMES + 1));
	auto zeros = std::string(max_digits - num.size(), '0');

	ret += zeros + num;

	return ret;
}

int main(int argc, char *argv[])
{
	static const struct option long_options[] = {
		{"version", no_argument, NULL, 0},
		{"help", no_argument, NULL, 0},
		{"gene-length", required_argument, NULL, 'l'},
		{"num-genes", required_argument, NULL, 'm'},
		{"num-genomes", required_argument, NULL, 'n'},
		{"out-dir", required_argument, NULL, 'o'},
		{"probability", required_argument, NULL, 'p'},
		{"seed", required_argument, NULL, 's'},
		{"verbose", no_argument, NULL, 'v'} // no comment
	};

	while (true) {
		int option_index = 0;
		int c = getopt_long(argc, argv, "e:l:m:n:o:p:s:t:v", long_options,
							&option_index);

		if (c == -1) {
			break;
		}

		try {
			switch (c) {
				case 0: {
					auto option_string =
						std::string(long_options[option_index].name);
					if (option_string == "help") {
						usage(EXIT_SUCCESS);
					}
					if (option_string == "version") {
						version();
					}
					return 0;
				}
				case 'e': {
					EVO_DISTANCE = std::stod(optarg);
					if (EVO_DISTANCE < 0.0 || EVO_DISTANCE > 1.0) {
						throw std::out_of_range(
							"evolutionary distance not in range [0,1]");
					}
					break;
				}
				case 'l': {
					GENE_LENGTH = std::stoul(optarg);
					break;
				}
				case 'n': {
					NUM_GENOMES = std::stoul(optarg);
					break;
				}
				case 'm': {
					NUM_GENES = std::stoul(optarg);
					break;
				}
				case 'o': {
					OUT_DIR = std::string(optarg) + "/";
					break;
				}
				case 'p': {
					PROBABILITY = std::stod(optarg);
					if (PROBABILITY < 0.0 || PROBABILITY > 1.0) {
						throw std::out_of_range(
							"probability not in range [0,1]");
					}
					break;
				}
				case 's': {
					SEED = std::stoul(optarg);
					break;
				}
				case 't': {
					THRESHOLD = std::stod(optarg);
					if (THRESHOLD < 0.0 || THRESHOLD > 1.0) {
						throw std::out_of_range(
							"core genome size not in range [0,1]");
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
		} catch (const std::exception &e) {
			warnx("error parsing -%c argument", c);
			throw;
		}
	}

	// setup RNG
	RNG = std::default_random_engine(SEED);
	auto rng_float_help = std::uniform_real_distribution<double>(0, 1);
	auto rng_float = bind(rng_float_help, RNG);

	// setup distribution of pan genome
	auto distribution = create_distribution();
	auto genome_numbers = std::vector<size_t>(NUM_GENOMES);
	std::iota(genome_numbers.begin(), genome_numbers.end(), 0);

	if (VERBOSE) {
		size_t i = 0;
		for (auto d : distribution) {
			std::cout << i++ << "\t" << d << "\n";
		}
	}

	// simulate pan genomes
	auto pan_genes = std::vector<pan_gene>();
	for (size_t ng = 0; ng < NUM_GENES; ng++) {
		auto random_float = rng_float();
		if (VERBOSE) {
			std::cout << random_float << "\n";
		}
		auto it = std::lower_bound(distribution.begin(), distribution.end(),
								   random_float);

		// FIXME: clean up next line
		auto pan_gene_size =
			NUM_GENOMES + 1 - std::distance(distribution.begin(), it);
		if (VERBOSE) {
			std::cout << pan_gene_size << "\n";
		}

		std::shuffle(genome_numbers.begin(), genome_numbers.end(), RNG);

		auto sub_set = std::vector<size_t>(
			genome_numbers.begin(), genome_numbers.begin() + pan_gene_size);

		// simulate genes
		pan_genes.emplace_back(sub_set);
	}

	// create output directory
	auto check =
		mkdir(OUT_DIR.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if (check && errno != EEXIST) {
		err(errno, "%s", OUT_DIR.c_str());
	}

	auto check_io = [](const std::ofstream &o, const std::string &n) {
		if (!o) {
			err(errno, "%s: couldn't write to file", n.c_str());
		}
	};

	{ // print reproducible information
		auto rep_file_name = OUT_DIR + "seed";
		auto rep_file = std::ofstream(rep_file_name);
		check_io(rep_file, rep_file_name);

		rep_file << "pangenomesim " VERSION << "\n";
		rep_file << SEED << "\n";
		// … arguments
		check_io(rep_file, rep_file_name);
	}

	// produce plenty of output
	std::sort(genome_numbers.begin(), genome_numbers.end());

	std::cout << "number of genomes:\t" << genome_numbers.size() << "\n";
	std::cout << "number of genes:\t" << pan_genes.size() << "\n";
	for (const auto &pg : pan_genes) {
		std::cout << "pg size:\t" << pg.size() << "\n";
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

		size_t chromosome = 0;
		for (const auto &pg : pan_genes) {
			ref_fasta_file << ">ref.chr" << ++chromosome << "\n"
						   << pg.reference() << std::endl;
			check_io(ref_fasta_file, ref_file_name);

			auto is_core = THRESHOLD * double(NUM_GENOMES) <= double(pg.size());
			auto &file = is_core ? cor_fasta_file : acc_fasta_file;
			auto &name = is_core ? cor_file_name : acc_file_name;

			file << ">ref.chr" << ++chromosome << "\n"
				 << pg.reference() << std::endl;
			check_io(file, name);
		}
	}

	// print derived genomes
	for (auto gn : genome_numbers) {
		auto file_name = OUT_DIR + genome_name(gn + 1) + ".fasta";
		auto fasta_file = std::ofstream(file_name);
		check_io(fasta_file, file_name);

		size_t chromosome = 0;
		for (const auto &pg : pan_genes) {
			auto idx = pg.find_genome_number(gn);
			if (idx == pg.size()) {
				continue; // not found
			}

			fasta_file << ">" << genome_name(gn + 1) << ".chr" << ++chromosome
					   << "\n"
					   << pg.derived(idx) << std::endl;
			check_io(fasta_file, file_name);
		}
	}

	{ // print MAF alignment
		auto reference_counter = 1;
		auto chromosome_counters = std::vector<size_t>(NUM_GENOMES, 1);
		auto maf_file_name = OUT_DIR + "alignment.maf";
		auto maf_file = std::ofstream(maf_file_name);
		check_io(maf_file, maf_file_name);

		maf_file << "##maf version=1\n";
		check_io(maf_file, maf_file_name);

		for (const auto &pg : pan_genes) {
			maf_file << "a\n";
			maf_file << "s ref.chr" << reference_counter++ << "\n"
					 << pg.reference() << std::endl;

			auto size = pg.size();
			for (size_t k = 0; k < size; k++) {
				auto gn = pg.genome_number(k);
				maf_file << "s " << genome_name(gn + 1) << ".chr"
						 << chromosome_counters[gn]++ << " " << pg.reference()
						 << std::endl;
			}

			maf_file << std::endl;
			check_io(maf_file, maf_file_name);
		}
	}

	return 0;
}

void usage(int exit_code)
{
	static const char str[] = {
		"usage: pangenomesim [OPTIONS...]\n"
		"Simulate a pan-genome.\n\n"
		"  -e, --evo-distance=FLOAT Set the evolutionary distance\n"
		"  -l, --gene-length=NUM    Set the gene length\n"
		"  -m, --num-genes=NUM      Set the number of genes\n"
		"  -n, --num-genomes=NUM    Set the number of genomes\n"
		"  -o, --out-dir=DIR        The directory to write files to\n"
		"  -p, --probability=FLOAT  Set simulation parameter\n"
		"  -s, --seed=NUM           Seed for the random number generator\n"
		"  -t, --threshold=FLOAT    Set the relative core genome size\n"
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
		"Copyright (C) 2017  Fabian Klötzl <kloetzl@evolbio.mpg.de>\n"
		"BSD License\n" //
	};
	printf("%s", str);
	exit(EXIT_SUCCESS);
}
