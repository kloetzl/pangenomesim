#pragma once
#include <string>
#include <vector>

#include "locus.h"

class evo_model
{
  protected:
  public:
	evo_model() = default;
	virtual ~evo_model() = default;

	virtual void parse_param(std::string, std::string);
	virtual std::string parameters() const;
	virtual void simulate() = 0;

	virtual std::vector<locus> get_reference() = 0;
	virtual std::vector<locus> get_core() = 0;
	virtual std::vector<locus> get_accessory() = 0;
	virtual std::vector<locus> get_genome(ssize_t) = 0;
	virtual std::vector<locus> get_locus(ssize_t) = 0;
};
