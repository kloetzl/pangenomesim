#pragma once
#include "evo_model.h"
#include <string>

class img_model : public evo_model
{
  public:
	img_model() = default;

	std::string parameters();
	void simulate();

	std::vector<locus> get_reference();
	std::vector<locus> get_core();
	std::vector<locus> get_accessory();
	std::vector<locus> get_genome(size_t);
	std::vector<locus> get_locus(size_t);
};
