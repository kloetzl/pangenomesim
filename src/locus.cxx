#include "locus.h"
#include <string>

locus::locus(std::string n, size_t gid, size_t lid)
	: nucl(n), genome_id(gid), locus_id(lid)
{
}

const std::string &locus::get_nucl() const
{
	return nucl;
}

size_t locus::get_genome_id() const
{
	return genome_id;
}

size_t locus::get_locus_id() const
{
	return locus_id;
}
