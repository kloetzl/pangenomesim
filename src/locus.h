#pragma once
#include "util.h"
#include <string>
#include <vector>

class locus
{
	std::string nucl = {};
	ssize_t genome_id = -1;
	ssize_t locus_id = -1;

  public:
	locus() = default;
	locus(std::string, ssize_t, ssize_t);
	locus(ssize_t, ssize_t, ssize_t);

	const std::string &get_nucl() const;
	ssize_t get_genome_id() const;
	ssize_t get_locus_id() const;
	void set_genome_id(ssize_t);

	locus mutate(double) const;
	static std::vector<locus> vector_mutate(const std::vector<locus> &, double);

	std::string to_fasta() const
	{
		auto ret = std::string(">");
		ret += genome_name(genome_id);
		ret += ".";
		ret += std::to_string(locus_id);
		ret += "\n";
		ret += nucl;
		ret += "\n";
		return ret;
	}
};
