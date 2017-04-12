#pragma once
#include <string>
#include <vector>

#include "locus.h"

class evo_model
{
  protected:
	size_t loci_length = 100;
	size_t num_loci = 5;
	size_t num_genomes = 3;
	size_t seed = 1729;

  public:
	evo_model() = default;
	virtual ~evo_model() = default;

	virtual void parse_param(std::string, std::string);
	virtual std::string parameters() const;
	virtual void simulate() = 0;

	virtual std::vector<locus> get_reference() = 0;
	virtual std::vector<locus> get_core() = 0;
	virtual std::vector<locus> get_accessory() = 0;
	virtual std::vector<locus> get_genome(size_t) = 0;
	virtual std::vector<locus> get_locus(size_t) = 0;
};
