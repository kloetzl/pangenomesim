#include "global.h"
#include "global.h"
#include "pangene.h"
#include "simulate.h"
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
#include <vector>

#ifndef VERSION
#define VERSION "v0.1"
#endif

void usage(int);
void version();

double EVO_DISTANCE = 0.01;
size_t GENE_LENGTH = 100;
size_t NUM_GENES = 5;
size_t NUM_GENOMES = 3;
std::string OUT_DIR = std::string("./");
double PROBABILITY = 0.8;
size_t SEED = 1729;
double THRESHOLD = 1.0;
std::default_random_engine RNG;
bool VERBOSE = false;

int main(int argc, char *argv[])
{
	static const struct option long_options[] = {
		{"version", no_argument, NULL, 0},
		{"help", no_argument, NULL, 0},
		{"evo-distance", required_argument, NULL, 'e'},
		{"gene-length", required_argument, NULL, 'l'},
		{"num-genes", required_argument, NULL, 'm'},
		{"num-genomes", required_argument, NULL, 'n'},
		{"out-dir", required_argument, NULL, 'o'},
		{"probability", required_argument, NULL, 'p'},
		{"seed", required_argument, NULL, 's'},
		{"threshold", required_argument, NULL, 't'},
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

	simulate();

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
		"Copyright (C) 2017,  Fabian Klötzl "
		"<fabian-pangenomesim@kloetzl.info>\n"
		"BSD 2-Clause License\n" //
	};
	printf("%s", str);
	exit(EXIT_SUCCESS);
}
