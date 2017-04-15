#pragma once
#include "evo_model.h"
#include <string>

class img_model : public evo_model
{
  protected:
	double img_theta = 0.1;
	double img_rho = 0.1;
	size_t img_core_size = 4;
	std::vector<locus> ref_core = {};
	std::vector<locus> ref_acc = {};
	std::vector<std::vector<locus>> loci = {};
	std::vector<std::vector<locus>> acc = {};

  public:
	img_model() = default;

	virtual void parse_param(std::string, std::string);
	std::string parameters() const;
	void simulate();

	std::vector<locus> get_reference();
	std::vector<locus> get_core();
	std::vector<locus> get_accessory();
	std::vector<locus> get_genome(size_t);
	std::vector<locus> get_locus(size_t);
};
