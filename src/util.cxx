#include "util.h"
#include "global.h"
#include <cmath>
#include <string>

extern "C" {
#include <err.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
}

std::string genome_name(size_t genome_number)
{
	auto ret = std::string("genome");
	auto num = std::to_string(genome_number);
	// TODO: get number of genomes?
	auto max_digits = std::ceil(std::log10(/*NUM_GENOMES*/ 100 + 1));
	auto zeros = std::string(max_digits - num.size(), '0');

	ret += zeros + num;

	return ret;
}

void mkpath(std::string path)
{
	const auto flags = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
	int check = ::mkdir(path.c_str(), flags);
	if (check == 0 || errno == EEXIST) {
		return; // success
	}
	if (errno != ENOENT) {
		err(errno, "%s", path.c_str());
	}

	auto parent = path.substr(0, path.find_last_of('/'));
	mkpath(parent);
	check = ::mkdir(path.c_str(), flags);
	if (check && errno != EEXIST) {
		err(errno, "%s", path.c_str());
	}
}