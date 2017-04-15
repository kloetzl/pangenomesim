#pragma once
#include <string>

std::string genome_name(ssize_t);
void mkpath(std::string);

size_t rand_int(size_t lower, size_t upper);
double rand_exp(double arg);
size_t rand_poisson(double arg);
