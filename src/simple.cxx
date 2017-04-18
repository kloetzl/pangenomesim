#include "simple.h"
#include "config.h"
#include "global.h"
#include "pangene.h"
#include "util.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

extern "C" {
#include <err.h>
#include <errno.h>
}

void simple_model::parse_param(std::string key, std::string value)
{

	if (key == "probability") {
		probability = std::stod(value);
		if (probability < 0.0 || probability > 1.0) {
			throw std::out_of_range("probability not in range [0,1]");
		}
	} else if (key == "threshold") {
		threshold = std::stod(value);
		if (threshold < 0.0 || threshold > 1.0) {
			throw std::out_of_range("core genome size not in range [0,1]");
		}
	} else if (key == "evo_distance") {
		evo_distance = std::stod(value);
		if (evo_distance < 0.0 || evo_distance > 1.0) {
			throw std::out_of_range("evolutionary distance not in range [0,1]");
		}
	} else {
		evo_model::parse_param(key, value);
	}
}

std::vector<double> simple_model::create_distribution()
{
	auto distribution = std::vector<double>(num_genomes + 1);

	auto tail = std::pow(1 - probability, num_genomes);
	auto factor = 1 / (1 - tail);
	std::cout << "factor\t" << factor << "\n";

	for (size_t i = 1; i < num_genomes + 1; i++) {
		distribution[i] =
			distribution[i - 1] +
			probability * std::pow(1 - probability, i - 1) * factor;
	}

	return distribution;
}

std::string simple_model::parameters() const
{
	auto str = std::stringstream();

	str << "pangenomesim " VERSION << "\n\n";
	str << "full options:\n";
	str << "--param evo-distance=" << evo_distance << "\n";
	str << "--param loci-length=" << loci_length << "\n";
	str << "--param num-loci=" << num_loci << "\n";
	str << "--param num-genomes=" << num_genomes << "\n";
	str << "--param probability=" << probability << "\n";
	str << "--param seed=" << seed << "\n";
	str << "--param threshold=" << threshold << "\n";

	return str.str();
}

std::vector<locus> simple_model::get_reference()
{
	auto ret = std::vector<locus>();

	size_t chromosome = 0;
	for (const auto &pg : pan_genes) {
		auto loc = locus(pg.reference(), -1, ++chromosome);
		ret.push_back(loc);
	}

	return ret;
}

std::vector<locus> simple_model::get_core()
{
	auto ret = std::vector<locus>();

	size_t chromosome = 0;
	for (const auto &pg : pan_genes) {
		++chromosome;

		if (threshold * double(num_genomes) <= double(pg.size())) {
			auto loc = locus(pg.reference(), -1, chromosome);
			ret.push_back(loc);
		}
	}

	return ret;
}

std::vector<locus> simple_model::get_accessory()
{
	auto ret = std::vector<locus>();

	size_t chromosome = 0;
	for (const auto &pg : pan_genes) {
		++chromosome;

		if (threshold * double(num_genomes) > double(pg.size())) {
			auto loc = locus(pg.reference(), -1, chromosome);
			ret.push_back(loc);
		}
	}

	return ret;
}

std::vector<locus> simple_model::get_genome(ssize_t gn)
{
	auto ret = std::vector<locus>();

	size_t chromosome = 0;
	for (const auto &pg : pan_genes) {
		auto idx = pg.find_genome_number(gn);
		if (idx == pg.size()) {
			continue; // not found
		}
		chromosome++;

		auto loc = locus(pg.derived(idx), gn, chromosome);
		ret.push_back(loc);
	}

	return ret;
}

std::vector<locus> simple_model::get_locus(ssize_t ln)
{
	auto ret = std::vector<locus>();
	const auto &pg = pan_genes.at(ln);

	auto loc = locus(pg.reference(), -1, ln);
	ret.push_back(loc);

	auto size = pg.size();
	for (size_t k = 0; k < size; k++) {
		auto gn = pg.genome_number(k);
		// todo: this is inefficient
		size_t chromosome = 0;
		for (size_t k = 0; k < ln; k++) {
			if (pan_genes[k].find_genome_number(gn) != pan_genes[k].size()) {
				chromosome++;
			}
		}
		loc = locus(pg.derived(k), gn, chromosome);
	}

	return ret;
}

void simple_model::simulate()
{
	// setup rng
	RNG = std::default_random_engine(seed);
	auto rng_float_help = std::uniform_real_distribution<double>(0, 1);
	auto rng_float = bind(rng_float_help, RNG);

	// setup distribution of pan genome
	auto distribution = create_distribution();
	auto genome_numbers = std::vector<size_t>(num_genomes);
	std::iota(genome_numbers.begin(), genome_numbers.end(), 0);

	if (VERBOSE) {
		size_t i = 0;
		for (auto d : distribution) {
			std::cout << i++ << "\t" << d << "\n";
		}
	}

	// simulate pan genomes
	pan_genes = std::vector<pan_gene>();
	for (size_t ng = 0; ng < num_loci; ng++) {
		auto random_float = rng_float();
		if (VERBOSE) {
			std::cout << random_float << "\n";
		}
		auto it = std::lower_bound(distribution.begin(), distribution.end(),
								   random_float);

		// fixme: clean up next line
		auto pan_gene_size =
			num_genomes + 1 - std::distance(distribution.begin(), it);
		if (VERBOSE) {
			std::cout << pan_gene_size << "\n";
		}

		std::shuffle(genome_numbers.begin(), genome_numbers.end(), RNG);

		auto sub_set = std::vector<size_t>(
			genome_numbers.begin(), genome_numbers.begin() + pan_gene_size);

		// simulate genes
		pan_genes.emplace_back(this, sub_set);
	}

	std::sort(genome_numbers.begin(), genome_numbers.end());

	if (VERBOSE) {
		std::cout << "number of genomes:\t" << genome_numbers.size() << "\n";
		std::cout << "number of genes:\t" << pan_genes.size() << "\n";
		for (const auto &pg : pan_genes) {
			std::cout << "pg size:\t" << pg.size() << "\n";
		}
	}

	auto genome_sizes = std::vector<size_t>(num_genomes);
}
