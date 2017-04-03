#pragma once

#include <string>
#include <vector>

class pan_gene
{
	std::string reference_sequence = "";
	std::vector<std::string> derived_sequences = {};
	std::vector<size_t> genome_numbers = {};

  public:
	pan_gene() = default;
	pan_gene(std::vector<size_t> genome_numbers);

	size_t size() const noexcept;
	const std::string &reference() const noexcept;
	const std::string &derived(size_t index) const;
	size_t genome_number(size_t index) const;
	size_t find_genome_number(size_t gn) const;
};
