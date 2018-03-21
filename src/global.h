#pragma once
#include "gene.h"
#include <random>
#include <string>

extern std::string OUT_DIR;
extern std::default_random_engine RNG;
extern bool VERBOSE;

extern gene root_gene(ssize_t, ssize_t);
