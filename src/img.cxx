#include "img.h"
#include "config.h"
#include "global.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

extern "C" {
#include <err.h>
}

void img_model::parse_param(std::string key, std::string value)
{
	if (key == "theta") {
		this->img_theta = std::stof(value);
		return;
	}
	if (key == "rho") {
		this->img_rho = std::stof(value);
		return;
	}
	if (key == "core-size") {
		this->img_core_size = std::stoul(value);
		return;
	}
	if (key == "gene-length") {
		this->gene_length = std::stoul(value);
		return;
	}
	if (key == "num-genomes") {
		this->num_genomes = std::stoul(value);
		return;
	}
	if (key == "seed") {
		this->seed = std::stoul(value);
		RNG = std::default_random_engine(seed);
		return;
	}
	if (key == "mut-rate") {
		this->mut_rate = std::stof(value);
		return;
	}
	errx(1, "unkown key %s.", key.c_str());
}

std::string img_model::parameters() const
{
	auto str = std::stringstream();

	str << "--param core-size=" << img_core_size << "\n";
	str << "--param gene-length=" << gene_length << "\n";
	str << "--param mut-rate=" << mut_rate << "\n";
	str << "--param num-genomes=" << num_genomes << "\n";
	str << "--param rho=" << img_rho << "\n";
	str << "--param seed=" << seed << "\n";
	str << "--param theta=" << img_theta << "\n";

	return str.str();
}

/**
 * @brief Simulate a coalescent.
 * @param n - The sample size.
 * @returns The nodes of the coalescent, allocated in a pool. Node 0 to n are
 * leaves, the rest are internal branches with the last being the root.
 */
std::vector<tree_node> create_coalescent(size_t n)
{
	auto pool = std::vector<tree_node>(2 * n - 1);

	generate_i(pool.begin(), n, [](size_t i) {
		return tree_node::leaf(i); //
	});

	// reorder the `indirect` references, rather than the pool itself.
	auto indirect = std::vector<tree_node *>();
	indirect.reserve(2 * n - 1);
	for (auto &tn : pool) {
		indirect.push_back(&tn);
	}

	// simulate topology
	for (size_t i = 0; i < n - 1; i++) {
		auto p = n + i;
		auto c = rand_int(0, n - i);
		indirect.at(p)->set_left(indirect[c]); // also sets [c].parent
		indirect[c] = indirect[n - i - 1];

		auto d = rand_int(0, n - i - 1);
		indirect[p]->set_right(indirect[d]);
		indirect[d] = indirect[p];
	}

	// add times
	auto t = 0.0;

	// simulate absolute times
	auto sample_size = n;
	for (size_t i = 0; i < n - 1; i++, sample_size--) {
		auto ss1 = sample_size * (sample_size - 1) / 2; // N_e ??
		t += rand_exp(ss1);
		pool[n + i].set_abs_time(t);
	}

	auto root_it = pool.end() - 1;
	// compute relative times
	std::for_each(pool.begin(), root_it, [](tree_node &self) {
		self.compute_rel_time(); //
	});

	return pool;
}

/**
 * @brief Compute distance matrix
 */
std::vector<std::vector<double>>
create_distmatrix(const std::vector<tree_node> &pool, double mut_rate)
{
	auto INF = std::numeric_limits<double>::infinity();
	auto m = pool.size();
	auto n = (pool.size() + 1) / 2;
	auto big = std::vector<double>(m * m, INF);
	auto B = [&](size_t x, size_t y) -> double & {
		return big[x * m + y]; // hack
	};

	for (size_t i = 0; i < m; i++) {
		B(i, i) = 0.0;
	}

	auto start_ptr = &pool[0];
	auto &root = top(pool);
	root.traverse([&](const tree_node &self) {
		if (self.has_parent()) {
			auto x = &self - start_ptr;
			auto y = &self.get_parent() - start_ptr;
			B(x, y) = B(y, x) = self.get_time() * mut_rate;
		}
	});

	// Floyd Warshall Algorithm
	for (size_t k = 0; k < m; k++) {
		for (size_t x = 0; x < m; x++) {
			for (size_t y = 0; y < m; y++) {
				auto tmp = std::min(B(x, y), B(x, k) + B(k, y));
				B(x, y) = B(y, x) = tmp;
			}
		}
	}

	auto base = std::vector<double>(n, 0.0);
	auto ret = std::vector<std::vector<double>>(n, base);
	for (size_t x = 0; x < n; x++) {
		for (size_t y = 0; y < n; y++) {
			ret[x][y] = ret[y][x] = B(x, y);
		}
	}

	return ret;
}

/**
 * @brief Simulate sequence evolution from the root.
 * @returns a vector of sequences sorted by genome ID.
 */
std::vector<gene> img_model::seq_from_root(const tree_node &root,
										   size_t gene_id)
{
	// create core sequences
	auto sample_size = num_genomes;
	auto leaves = std::vector<gene>(sample_size);
	auto stack = std::vector<gene>();
	stack.reserve(sample_size);
	auto seq = gene(gene_length, -1, gene_id); // random root seq
	stack.push_back(seq);
	ref_core.push_back(seq);
	root.traverse(
		[&stack, this](const tree_node &self) {
			if (self.has_parent()) {
				auto seq = top(stack).mutate(self.get_time() * mut_rate);
				stack.push_back(seq);
			}
		},
		[&stack, &leaves](const tree_node &self) {
			if (self.is_leaf()) {
				// set genome id
				top(stack).set_genome_id(self.get_index());
				leaves[self.get_index()] = top(stack);
			}
		},
		[&stack](const tree_node &) {
			stack.pop_back(); //
		});

	return leaves;
}

/**
 * @brief do the work.
 */
void img_model::simulate()
{
	if (seed == 0) {
		std::random_device rd;
		seed = rd();
		RNG = std::default_random_engine(seed);
	}

	// generate coalescent
	coalescent = create_coalescent(num_genomes);
	distmatrix = create_distmatrix(coalescent, mut_rate);
	auto &root = top(coalescent);

	// create core sequences
	generate_i(
		std::back_inserter(cor_genes), img_core_size,
		[this, &root](size_t gene_id) { return seq_from_root(root, gene_id); });

	// generate pan genome and create sequences
	auto acc_genes = std::vector<std::vector<gene>>();

	auto start = std::vector<gene>();
	auto start_size = rand_poisson(img_theta / img_rho);
	start.reserve(start_size);
	acc_genes.resize(start_size);
	generate_i(std::back_inserter(start), start_size, [=](size_t gene_id) {
		return gene(gene_length, -1, gene_id + img_core_size);
	});

	auto gene_id = start_size + img_core_size;
	ref_acc = start; // copy

	auto stack = std::vector<std::vector<gene>>();
	stack.reserve(num_genomes);
	stack.push_back(start);

	root.traverse(
		[&stack, &that = *this, &acc_genes, &gene_id ](const tree_node &self) {
			if (!self.has_parent()) {
				return; // root
			}

			// simulate evolution
			auto neu = gene::vector_mutate(top(stack),
										   that.mut_rate * self.get_time());
			auto time = 0.0;
			while (time < self.get_time()) {
				auto time_to_go = self.get_time() - time;
				auto time_to_gain = rand_exp(that.img_theta);
				auto time_to_loss = rand_exp(that.img_rho * neu.size());
				if (time_to_gain > time_to_go && time_to_loss > time_to_go) {
					break;
				}
				// gain or loss event
				/* This works because gene-gain/loss is a Poisson process
				 * and as such the individual events are exponentially
				 * distributed. This means the distribution is memoryless and
				 * thus past events can be ignored.
				 */
				if (time_to_gain < time_to_loss || neu.empty()) {
					neu.emplace_back(that.gene_length, -1, gene_id++);
					acc_genes.resize(acc_genes.size() + 1);
					that.ref_acc.push_back(top(neu));
					time += time_to_gain;
				} else {
					assert(neu.size() != 0);
					auto loser = rand_int(0, neu.size());
					// pick one
					std::swap(neu[loser], top(neu));
					neu.pop_back();
					time += time_to_loss;
				}
			}

			stack.push_back(neu);
		},
		[&stack, &acc_genes,
		 img_core_size = this->img_core_size ](const tree_node &self) {
			if (self.is_branch()) {
				return;
			}

			auto genome_id = self.get_index();
			for (auto &loc : top(stack)) {
				loc.set_genome_id(genome_id);
				acc_genes[loc.get_gene_id() - img_core_size].push_back(loc);
			}
		},
		[&stack](const tree_node &) {
			stack.pop_back(); //
		});

	assert(stack.empty());

	// filter out empty genes which were fully lost
	assert(ref_acc.size() == acc_genes.size());
	auto it =
		remove_if_i(ref_acc.begin(), ref_acc.end(),
					[&acc_genes](size_t i) { return acc_genes.at(i).empty(); });
	ref_acc.erase(it, ref_acc.end());

	auto is_empty = [](const std::vector<gene> &set) { return set.empty(); };

	// erase-remove-idiom
	auto empty_it =
		std::remove_if(acc_genes.begin(), acc_genes.end(), is_empty);
	acc_genes.erase(empty_it, acc_genes.end());

	// move genes that were not lost on any lineages from accessory to core.
	auto is_core = [=](const std::vector<gene> &set) {
		return set.size() == num_genomes;
	};
	for (size_t i = 0; i < acc_genes.size(); i++) {
		if (is_core(acc_genes[i])) {
			cor_genes.push_back(acc_genes[i]);
			ref_core.push_back(ref_acc[i]);
		}
	}

	auto core_it = std::remove_if(acc_genes.begin(), acc_genes.end(), is_core);
	acc_genes.erase(core_it, acc_genes.end());

	this->acc_genes = acc_genes;
}

std::vector<gene> img_model::get_reference()
{
	auto ret = ref_core;
	ret.reserve(ref_core.size() + ref_acc.size());

	std::copy(ref_acc.begin(), ref_acc.end(), std::back_inserter(ret));

	return ret;
}

std::vector<gene> img_model::get_core()
{
	return ref_core; // implicit copy
}

std::vector<gene> img_model::get_accessory()
{
	return ref_acc; // implicit copy
}

std::vector<gene> img_model::get_genome(ssize_t genome_id)
{
	auto ret = std::vector<gene>();
	auto inserter = std::back_inserter(ret);

	auto gid_filter = [&genome_id](const auto &loc) {
		return loc.get_genome_id() == genome_id;
	};

	for (const auto &gene_set : cor_genes) {
		std::copy_if(gene_set.begin(), gene_set.end(), inserter, gid_filter);
	}

	for (const auto &gene_set : acc_genes) {
		std::copy_if(gene_set.begin(), gene_set.end(), inserter, gid_filter);
	}

	return ret;
}

std::vector<gene> img_model::get_gene(ssize_t some_number)
{
	if (some_number < cor_genes.size()) {
		return cor_genes.at(some_number); // throws on some_number < 0
	} else {
		return acc_genes.at(some_number - cor_genes.size()); // may throw
	}
}
