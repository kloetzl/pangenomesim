#include "evo_model.h"
#include "global.h"
// #include "img.h"
#include "locus.h"
#include "simple.h"
#include "util.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <functional>
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
size_t SEED = 1729;
std::default_random_engine RNG;
bool VERBOSE = false;

int main(int argc, char *argv[])
{
	static const struct option long_options[] = {
		{"version", no_argument, NULL, 0},
		{"help", no_argument, NULL, 0},
		{"param", required_argument, NULL, 0},
		{"model", required_argument, NULL, 0},
		{"out-dir", required_argument, NULL, 'o'},
		{"verbose", no_argument, NULL, 'v'},
		{0, 0, 0, 0} // no comment
	};

	// polymorphism, yay!
	evo_model *model = nullptr;

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
				} else if (option_string == "param") {
					if (!model) {
						errx(1, "set an evolutionary model via --model first");
					}

					auto arg_string = std::string(optarg);

					// implement `getsubsopt(3)` like functionality
					auto pattern = std::string("^(\\w+)(?:=(\\w+))?,?");
					auto r = std::regex(pattern);

					auto m = std::smatch();
					while (std::regex_search(arg_string, m, r)) {
						auto key = m[1];
						auto value = m.size() == 3 ? m[2] : std::string();
						std::cerr << key << ":" << value << std::endl;

						// may throw
						model->parse_param(key, value);
						arg_string = m.suffix();
					}

					if (!arg_string.empty()) {
						errx(1, "invalid parameter %s", optarg);
					}

					break;
				} else if (option_string == "model") {
					auto model_name = std::string(optarg);
					if (model_name == "simple") {
						model = new simple_model();
					// } else if (model_name == "IMG") {
						// model = new img_model();
					} else {
						errx(1, "unknown model: %s", model_name.c_str());
					}
				} else {
					return 1;
				}
				break;
			}
			case 'o': {
				OUT_DIR = std::string(optarg) + "/";
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

	if (!model) {
		errx(1, "no evolutionary model selected");
	}
	model->simulate(); // do the work

	auto check_io = [](const std::ofstream &o, const std::string &n) {
		if (!o) {
			err(errno, "%s: couldn't write to file", n.c_str());
		}
	};

	// TODO: output
	auto rep_file_name = OUT_DIR + "reproducible.seed";
	auto rep_file = std::ofstream(rep_file_name);
	check_io(rep_file, rep_file_name);

	rep_file << "pangenomesim " VERSION << "\n\n";
	rep_file << "full options:\n";
	rep_file << model->parameters() << std::endl;

	check_io(rep_file, rep_file_name);

	delete model;
	return 0;
}

void usage(int exit_code)
{
	static const char str[] = {
		"usage: pangenomesim [OPTIONS...]\n"
		"Simulate a pan-genome.\n\n"
		"  -o, --out-dir=DIR        The directory to write files to\n"
		"  -p, --probability=FLOAT  Set simulation parameter\n"
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
