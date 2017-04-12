#include "evo_model.h"
#include <exception>
#include <stdexcept>
#include <string>
#include <sstream>

void evo_model::parse_param(std::string key, std::string value)
{
	if (key == "loci_length") {
		this->loci_length = std::stoul(value);
		return;
	}
	if (key == "num_loci") {
		this->num_loci = std::stoul(value);
		return;
	}
	if (key == "num_genomes") {
		this->num_genomes = std::stoul(value);
		return;
	}
	if (key == "seed") {
		this->seed = std::stoul(value);
		return;
	}
	throw std::invalid_argument(std::string("unkown parameter ") + key);
}


std::string evo_model::parameters() const
{
	auto str = std::stringstream();

	str << "--param loci_length=" << loci_length << "\n";
	str << "--param num_loci=" << num_loci << "\n";
	str << "--param num_genomes=" << num_genomes << "\n";
	str << "--param seed=" << seed << "\n";

	return str.str();
}
