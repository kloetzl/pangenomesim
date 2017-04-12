#pragma once
#include "evo_model.h"
#include "locus.h"
#include "pangene.h"
#include <string>
#include <vector>

class simple_model : public evo_model
{
  protected:
	double evo_distance = 0.01;
	double probability = 0.8;
	double threshold = 1.0;
	std::vector<pan_gene> pan_genes = {};

  public:
	simple_model() = default;

	std::string parameters();
	virtual void parse_param(std::string, std::string);
	void simulate();

	std::vector<locus> get_reference();
	std::vector<locus> get_core();
	std::vector<locus> get_accessory();
	std::vector<locus> get_genome(size_t);
	std::vector<locus> get_locus(size_t);

  private:
	std::vector<double> create_distribution();

	friend class pan_gene;
};
