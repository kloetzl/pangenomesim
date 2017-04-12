#include "pangene.h"
#include "global.h"
#include "simple.h"
#include <algorithm>
#include <functional>
#include <string>
#include <vector>

static auto ACGT = "ACGT";
static auto NO_A = "CGT";
static auto NO_C = "AGT";
static auto NO_G = "ACT";
static auto NO_T = "ACG";

std::string pan_gene::random_seq() const
{
	auto ret = std::string();
	ret.reserve(m->loci_length);

	auto random_base_dist = std::uniform_int_distribution<int>(0, 3);
	auto random_base = [&] { return ACGT[random_base_dist(RNG)]; };

	for (size_t i = 0; i < m->loci_length; i++) {
		ret += random_base();
	}

	return ret;
}

std::string pan_gene::derive_seq(const std::string &reference) const
{
	auto ret = std::string(reference);

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

	double nucleotides = m->loci_length;
	double mutations = m->evo_distance * nucleotides;
	for (auto &c : ret) {
		if (mut() < mutations / nucleotides) {
			c = mutate(c);
			mutations--;
		}
		nucleotides--;
	}

	return ret;
}

pan_gene::pan_gene(const simple_model *mm, std::vector<size_t> __genome_numbers)
	: m(mm), genome_numbers(__genome_numbers)
{
	reference_sequence = random_seq();
	for (auto _ : genome_numbers) {
		derived_sequences.push_back(derive_seq(reference_sequence));
	}
}

size_t pan_gene::size() const noexcept
{
	return genome_numbers.size();
}

auto pan_gene::reference() const noexcept -> const std::string &
{
	return reference_sequence;
}

auto pan_gene::derived(size_t index) const -> const std::string &
{
	return derived_sequences.at(index);
}

auto pan_gene::genome_number(size_t index) const -> size_t
{
	return genome_numbers.at(index);
}

auto pan_gene::find_genome_number(size_t gn) const -> size_t
{
	auto it = std::find(genome_numbers.begin(), genome_numbers.end(), gn);
	return std::distance(genome_numbers.begin(), it);
}
