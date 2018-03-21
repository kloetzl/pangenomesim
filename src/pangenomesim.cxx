#include "config.h"
#include "gene.h"
#include "global.h"
#include "img.h"
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
#include "pfasta.h"
#include <err.h>
#include <errno.h>
#include <fcntl.h>
#include <getopt.h>
#include <sys/stat.h>
#include <unistd.h>
}

void usage(int);
void version();

std::string OUT_DIR = std::string("./");
std::default_random_engine RNG;
bool VERBOSE = false;
bool reference = false;
std::vector<gene> reference_pool{};

std::vector<gene> read_pool(std::string file_name);

int main(int argc, char *argv[])
{
	static const struct option long_options[] = {
		{"version", no_argument, NULL, 0},
		{"help", no_argument, NULL, 0},
		{"param", required_argument, NULL, 'p'},
		{"reference", required_argument, NULL, 'r'},
		{"out-dir", required_argument, NULL, 'o'},
		{"verbose", no_argument, NULL, 'v'},
		{0, 0, 0, 0} // no comment
	};

	auto model = img_model();

	if (argc == 1) usage(EXIT_SUCCESS);

	while (true) {
		int option_index = 0;
		int c = getopt_long(argc, argv, "o:p:r:v", long_options, &option_index);

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
				auto pattern = std::string("^([\\w\\-]+)(?:=([a-z_0-9.]+))?,?");
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
			case 'r': {
				reference = true;
				reference_pool = read_pool(optarg);
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

		auto genes = model.get_genome(gn);
		for (const auto &loc : genes) {
			fasta_file << loc.to_fasta();
		}

		check_io(fasta_file, file_name);
	}

	{ // compute gene frequency spectrum
		auto num_genes = model.get_num_genes();
		auto gfs = std::vector<size_t>(model.get_num_genomes() + 1);
		for (size_t i = 0; i < num_genes; i++) {
			auto that = model.get_gene(i);
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

		auto gene_length = model.get_gene_length();

		auto num_genomes = model.get_num_genomes();
		// compute sequence sizes
		auto seq_sizes = std::vector<size_t>(num_genomes);
		for (size_t i = 0; i < num_genomes; i++) {
			seq_sizes[i] = model.get_genome(i).size() * gene_length;
		}

		auto ref = model.get_reference();
		auto ref_size = ref.size() * gene_length;

		auto s_line = [&gene_length, &maf_file, &seq_sizes](const auto &loc) {
			auto gid = loc.get_genome_id();
			maf_file << "s " << genome_name(gid) // genome ID
					 << "." << loc.get_gene_id() // contig ID
					 << " 0"					 // start
					 << " " << gene_length		 // size
					 << " +"					 // strand
					 << " "
					 << seq_sizes[gid] // size of the entire source sequence
					 << " " << loc.get_nucl() // sequence
					 << "\n";
		};

		for (size_t i = 0; i < model.get_num_genes(); i++) {
			maf_file << "a\n";
			// print reference
			maf_file << "s " << genome_name(-1)				// genome ID
					 << "." << ref[i].get_gene_id() << " 0" // start
					 << " " << gene_length					// size
					 << " +"								// strand
					 << " " << ref_size						// size s.a.
					 << " " << ref[i].get_nucl()			// sequence
					 << "\n";

			auto tmp = model.get_gene(i);
			std::for_each(tmp.begin(), tmp.end(), s_line);

			maf_file << std::endl;
		}

		check_io(maf_file, file_name);
	}

	return 0;
}

std::vector<gene> read_pool(std::string s_file_name)
{
	auto ret = std::vector<gene>();
	const char *file_name = s_file_name.c_str();

	int file_descriptor =
		s_file_name != "-" ? open(file_name, O_RDONLY) : STDIN_FILENO;

	if (file_descriptor < 0) {
		err(errno, "%s", file_name);
	}

	int l;
	pfasta_file pf;

	if ((l = pfasta_parse(&pf, file_descriptor)) != 0) {
		errx(1, "%s: %s", file_name, pfasta_strerror(&pf));
	}

	int counter = 0;
	pfasta_seq ps;
	while ((l = pfasta_read(&pf, &ps)) == 0) {
		// convert to gene
		ret.emplace_back(ps.seq, ps.name, -1, counter++);
		pfasta_seq_free(&ps);
	}

	if (l < 0) {
		errx(1, "%s: %s", file_name, pfasta_strerror(&pf));
	}

	pfasta_free(&pf);
	close(file_descriptor);

	std::reverse(ret.begin(), ret.end()); // hack
	return ret;
}

gene root_gene(ssize_t gene_length, ssize_t gene_id)
{
	if (!reference) return gene(gene_length, -1, gene_id);

	if (reference_pool.empty()) {
		reference = false;
		std::cerr << "reference gene pool exhausted, switching to random\n";
		return gene(gene_length, -1, gene_id);
	}

	auto ret = reference_pool.back();
	reference_pool.pop_back();
	return ret;
}

void usage(int exit_code)
{
	static const char str[] = {
		"usage: pangenomesim [OPTIONS...]\n"
		"Simulate a pangenome according to the 'infinitely many genes' "
		"model.\n\n"
		"Options:\n"
		"  -o, --out-dir DIR        The directory to write files to\n"
		"  -p, --param KEY=VALUE    Set simulation parameter\n"
		"  -r, --reference FASTA    Set reference genes\n"
		"  -v, --verbose            Print additional information\n"
		"      --help               Display help and exit\n"
		"      --version            Output version information and exit\n\n"
		"Simulation Parameters:\n"
		"  core-size=INT            Minimum size of the core genome\n"
		"  gene-length=INT          Length of a gene\n"
		"  mut-rate=FLOAT           Substitution rate\n"
		"  num-genomes=INT          Number of genomes\n"
		"  rho=FLOAT                Rate of gene loss\n"
		"  seed=INT                 Seed for random number generator (use 0 "
		"for random)\n"
		"  theta=FLOAT              Rate of gene gain\n" //
	};

	fprintf(exit_code == EXIT_SUCCESS ? stdout : stderr, "%s", str);
	exit(exit_code);
}

void version()
{
	static const char str[] = {
		"pangenomesim " VERSION "\n"
		"Copyright (C) 2017 - 2018,  Fabian Klötzl "
		"<fabian-pangenomesim@kloetzl.info>\n"
		"BSD 2-Clause License\n" //
	};
	printf("%s", str);
	exit(EXIT_SUCCESS);
}
