#include "locus.h"
#include "global.h"
#include <functional>
#include <random>
#include <string>

locus::locus(std::string n, ssize_t gid, ssize_t lid)
	: nucl(n), genome_id(gid), locus_id(lid)
{
}

locus::locus(ssize_t length, ssize_t gid, ssize_t lid)
	: nucl(""), genome_id(gid), locus_id(lid)
{
	nucl.reserve(length);

	static const auto ACGT = "ACGT";
	auto random_base_dist = std::uniform_int_distribution<int>(0, 3);
	auto random_base = [&] { return ACGT[random_base_dist(RNG)]; };

	for (ssize_t i = 0; i < length; i++) {
		nucl += random_base();
	}
}

locus locus::mutate(double rate) const
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

const std::string &locus::get_nucl() const
{
	return nucl;
}

ssize_t locus::get_genome_id() const
{
	return genome_id;
}

ssize_t locus::get_locus_id() const
{
	return locus_id;
}

void locus::set_genome_id(ssize_t index)
{
	genome_id = index;
}
