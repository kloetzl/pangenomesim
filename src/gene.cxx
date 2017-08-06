#include "gene.h"
#include "global.h"
#include <functional>
#include <random>
#include <string>
#include <vector>

gene::gene(std::string n, ssize_t gid, ssize_t lid)
	: nucl(n), genome_id(gid), gene_id(lid)
{
}

gene::gene(ssize_t length, ssize_t gid, ssize_t lid)
	: nucl(""), genome_id(gid), gene_id(lid)
{
	nucl.reserve(length);

	static const auto ACGT = "ACGT";
	auto random_base_dist = std::uniform_int_distribution<int>(0, 3);
	auto random_base = [&] { return ACGT[random_base_dist(RNG)]; };

	for (ssize_t i = 0; i < length; i++) {
		nucl += random_base();
	}
}

gene gene::mutate(double rate) const
{
	auto ret = *this; // copy

	static auto NO_A = "CGT";
	static auto NO_C = "AGT";
	static auto NO_G = "ACT";
	static auto NO_T = "ACG";

	auto mut_dist = std::uniform_real_distribution<double>(0, 1);
	auto mut = bind(mut_dist, RNG);
	auto mut_acgt = std::uniform_int_distribution<int>(0, 2);
	auto mutate = [&](char c) {
		int idx = mut_acgt(RNG);
		switch (c) {
			case 'A': return NO_A[idx];
			case 'C': return NO_C[idx];
			case 'G': return NO_G[idx];
			case 'T': return NO_T[idx];
			default: return 'X';
		}
	};

	double nucleotides = nucl.size();
	double mutations = rate * nucleotides;
	for (auto &c : ret.nucl) {
		if (mut() < mutations / nucleotides) {
			c = mutate(c);
			mutations--;
		}
		nucleotides--;
	}

	return ret;
}

std::vector<gene> gene::vector_mutate(const std::vector<gene> &set,
										double rate)
{
	auto ret = std::vector<gene>();
	ret.reserve(set.size());
	for (const auto &loc : set) {
		ret.push_back(loc.mutate(rate));
	}
	return ret;
}

const std::string &gene::get_nucl() const
{
	return nucl;
}

ssize_t gene::get_genome_id() const
{
	return genome_id;
}

ssize_t gene::get_gene_id() const
{
	return gene_id;
}

void gene::set_genome_id(ssize_t index)
{
	genome_id = index;
}
