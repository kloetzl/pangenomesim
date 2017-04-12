#pragma once

#include <string>
#include <vector>

// forward declaration
class simple_model;

class pan_gene
{
	const simple_model *m = nullptr;
	std::string reference_sequence = "";
	std::vector<std::string> derived_sequences = {};
	std::vector<size_t> genome_numbers = {};

  public:
	pan_gene() = default;
	pan_gene(const simple_model *, std::vector<size_t> genome_numbers);

	size_t size() const noexcept;
	const std::string &reference() const noexcept;
	const std::string &derived(size_t index) const;
	size_t genome_number(size_t index) const;
	size_t find_genome_number(size_t gn) const;

  private:
	std::string random_seq() const;
	std::string derive_seq(const std::string &) const;
};
