#pragma once
#include <string>

class locus
{
	std::string nucl = {};
	size_t genome_id = 0;
	size_t locus_id = 0;

  public:
	locus() = default;
	locus(std::string, size_t, size_t);

	const std::string &get_nucl() const;
	size_t get_genome_id() const;
	size_t get_locus_id() const;
};
