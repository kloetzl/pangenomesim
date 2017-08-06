#pragma once
#include "util.h"
#include <string>
#include <vector>

class gene
{
	std::string nucl = {};
	ssize_t genome_id = -1;
	ssize_t gene_id = -1;

  public:
	gene() = default;
	gene(std::string, ssize_t, ssize_t);
	gene(ssize_t, ssize_t, ssize_t);

	const std::string &get_nucl() const;
	ssize_t get_genome_id() const;
	ssize_t get_gene_id() const;
	void set_genome_id(ssize_t);

	gene mutate(double) const;
	static std::vector<gene> vector_mutate(const std::vector<gene> &, double);

	std::string to_fasta() const
	{
		auto ret = std::string(">");
		ret += genome_name(genome_id);
		ret += ".";
		ret += std::to_string(gene_id);
		ret += "\n";
		ret += nucl;
		ret += "\n";
		return ret;
	}
};
