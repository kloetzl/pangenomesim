#pragma once
#include <string>

std::string genome_name(ssize_t);
void mkpath(std::string);

size_t rand_int(size_t lower, size_t upper);
double rand_exp(double arg);
size_t rand_poisson(double arg);

template <typename OutputIt, typename Size, typename Func>
OutputIt generate_i(OutputIt first, Size count, Func g)
{
	for (Size i = 0; i < count; ++i) {
		*first = g(i);
		first++;
	}
	return first;
}
