#include "config.h"
#include "global.h"
#include "pangene.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <vector>

extern "C" {
#include <err.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
}

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

void mkpath(const std::string &path)
{
	const auto flags = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
	int check = ::mkdir(path.c_str(), flags);
	if (check == 0 || errno == EEXIST) {
		return; // success
	}
	if (errno != ENOENT) {
		err(errno, "%s", path.c_str());
	}
	auto parent = path.substr(0, path.find_last_of('/'));
	mkpath(parent);
	check = ::mkdir(path.c_str(), flags);
	if (check) {
		err(errno, "%s", path.c_str());
	}
}

void simulate()
{
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
		auto rep_file_name = OUT_DIR + "reproducible.seed";
		auto rep_file = std::ofstream(rep_file_name);
		check_io(rep_file, rep_file_name);

		rep_file << "pangenomesim " VERSION << "\n\n";
		rep_file << "full options:\n";
		rep_file << "--evo-distance=" << EVO_DISTANCE << "\n";
		rep_file << "--gene-length=" << GENE_LENGTH << "\n";
		rep_file << "--num-genes=" << NUM_GENES << "\n";
		rep_file << "--num-genomes=" << NUM_GENOMES << "\n";
		rep_file << "--out-dir=" << OUT_DIR << "\n";
		rep_file << "--probability=" << PROBABILITY << "\n";
		rep_file << "--seed=" << SEED << "\n";
		rep_file << "--threshold=" << THRESHOLD << "\n";

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

	auto genome_sizes = std::vector<size_t>(NUM_GENOMES);

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

		genome_sizes[gn] = chromosome * GENE_LENGTH;
	}

	{ // print MAF alignment
		auto reference_counter = 0;
		auto chromosome_counters = std::vector<size_t>(NUM_GENOMES);
		auto maf_file_name = OUT_DIR + "alignment.maf";
		auto maf_file = std::ofstream(maf_file_name);
		check_io(maf_file, maf_file_name);

		maf_file << "##maf version=1\n";
		check_io(maf_file, maf_file_name);

		/* MAF format
		 *
		 * a
		 * s hg16.chr17    876234      13    +         6787  acgtagc
		 * s mm.chr12      872432      13    -         2345  TGCATAA
		 * s <name>.<chr>  <start> <size> <strand> <srcSize>  <seq>
		 */

		for (const auto &pg : pan_genes) {
			auto src_size = NUM_GENES * GENE_LENGTH;
			maf_file << "a\n";
			maf_file << "s ref.chr" << reference_counter << "\t" //
					 << reference_counter * GENE_LENGTH			 //
					 << "\t" << GENE_LENGTH						 //
					 << "\t" << src_size						 //
					 << "\t" << pg.reference() << std::endl;

			reference_counter++;

			auto size = pg.size();
			for (size_t k = 0; k < size; k++) {
				auto gn = pg.genome_number(k);
				maf_file << "s " << genome_name(gn + 1) << ".chr"
						 << (chromosome_counters[gn] + 1)				  //
						 << "\t" << chromosome_counters[gn] * GENE_LENGTH //
						 << "\t" << GENE_LENGTH							  //
						 << "\t" << genome_sizes[gn]					  //
						 << "\t" << pg.reference() << std::endl;

				chromosome_counters[gn]++;
			}

			maf_file << std::endl;
			check_io(maf_file, maf_file_name);
		}
	}
}