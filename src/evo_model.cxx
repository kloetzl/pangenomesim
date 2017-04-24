#include "evo_model.h"
#include "global.h"
#include <exception>
#include <sstream>
#include <stdexcept>
#include <string>

void evo_model::parse_param(std::string key, std::string value)
{

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
