#include "img.h"
#include "config.h"
#include "global.h"
#include <sstream>
#include <string>

void img_model::parse_param(std::string key, std::string value)
{
	if (key == "theta") {
		this->theta = std::stoul(value);
		return;
	}
	if (key == "rho") {
		this->rho = std::stoul(value);
		return;
	}
	evo_model::parse_param(key, value);
}

std::string img_model::parameters()
{
	auto str = std::stringstream();

	str << "pangenomesim " << VERSION << "\n\n";
	str << "full options:\n";
	str << "--param gene-length=" << loci_length << "\n";
	str << "--param num-genes=" << num_loci << "\n";
	str << "--param num-genomes=" << num_genomes << "\n";
	str << "--out-dir=" << OUT_DIR << "\n";
	str << "--param seed=" << seed << "\n";

	return str.str();
}

void img_model::simulate()
{
	// STUB
}

// STUB
std::vector<locus> img_model::get_reference()
{
	return std::vector<locus>();
}

// STUB
std::vector<locus> img_model::get_core()
{
	return std::vector<locus>();
}

// STUB
std::vector<locus> img_model::get_accessory()
{
	return std::vector<locus>();
}

// STUB
std::vector<locus> img_model::get_genome(size_t)
{
	return std::vector<locus>();
}

// STUB
std::vector<locus> img_model::get_locus(size_t)
{
	return std::vector<locus>();
}
