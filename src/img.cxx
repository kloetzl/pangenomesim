#include "img.h"
#include "config.h"
#include "global.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>

void img_model::parse_param(std::string key, std::string value)
{
	if (key == "img_theta") {
		this->img_theta = std::stof(value);
		return;
	}
	if (key == "img_rho") {
		this->img_rho = std::stof(value);
		return;
	}
	if (key == "img_core_size") {
		this->img_core_size = std::stoul(value);
		return;
	}
	if (key == "loci_length") {
		this->loci_length = std::stoul(value);
		return;
	}
	if (key == "num_genomes") {
		this->num_genomes = std::stoul(value);
		return;
	}
	if (key == "seed") {
		this->seed = std::stoul(value);
		RNG = std::default_random_engine(seed);
		return;
	}
}

std::string img_model::parameters() const
{
	auto str = std::stringstream();

	str << "pangenomesim " << VERSION << "\n\n";
	str << "full options:\n";
	str << "--param loci_length=" << loci_length << "\n";
	str << "--param num_genomes=" << num_genomes << "\n";
	str << "--param img_theta=" << img_theta << "\n";
	str << "--param img_rho=" << img_rho << "\n";
	str << "--param seed=" << seed << "\n";

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
		auto branch_length =
			self.get_parent().get_abs_time() - self.get_abs_time();
		self.set_time(branch_length);
	});

	return pool;
}

/**
 * @brief Simulate sequence evolution from the root.
 * @returns a vector of sequences sorted by genome ID.
 */
std::vector<locus> img_model::seq_from_root(const tree_node &root,
											size_t locus_id, double rate)
{
	// create core sequences
	auto sample_size = num_genomes;
	auto leaves = std::vector<locus>(sample_size);
	auto stack = std::vector<locus>();
	stack.reserve(sample_size);
	auto seq = locus(loci_length, -1, locus_id); // random root seq
	stack.push_back(seq);
	ref_core.push_back(seq);
	root.traverse(
		[&stack, &rate](const tree_node &self) {
			if (self.has_parent()) {
				auto seq = (stack.end() - 1)->mutate(self.get_time() * rate);
				stack.push_back(seq);
			}
		},
		[&stack, &leaves](const tree_node &self) {
			if (self.is_leaf()) {
				// set genome id
				auto &top = *(stack.end() - 1);
				top.set_genome_id(self.get_index());
				leaves[self.get_index()] = top;
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
	// generate coalescent
	coalescent = create_coalescent(num_genomes);
	auto &pool = coalescent;
	auto &root = *(pool.end() - 1);

	auto rate = 0.01;

	// create core sequences
	generate_i(std::back_inserter(cor_loci), img_core_size,
			   [&rate, this, &root](size_t locus_id) {
				   return seq_from_root(root, locus_id, rate);
			   });

	using loci_set = std::vector<locus>;

	// generate pan genome and create sequences
	auto acc_loci = std::vector<loci_set>();

	auto start = std::vector<locus>();
	auto start_size = rand_poisson(img_theta / img_rho);
	start.reserve(start_size);
	acc_loci.resize(start_size);
	generate_i(std::back_inserter(start), start_size, [=](size_t locus_id) {
		return locus(loci_length, -1, locus_id + img_core_size);
	});

	auto locus_id = start_size + img_core_size;
	ref_acc = start; // copy

	auto stack = std::vector<loci_set>();
	stack.reserve(num_genomes);
	stack.push_back(start);

	root.traverse(
		[&stack, &rate, &that = *this, &acc_loci,
		 &locus_id ](const tree_node &self) {
			if (!self.has_parent()) {
				return; // root
			}

			const auto &top = *(stack.end() - 1);
			// simulate evolution
			auto neu = locus::mutate_set(top, rate * self.get_time());
			auto time = 0.0;
			while (time < self.get_time()) {
				auto time_to_go = self.get_time() - time;
				auto time_to_gain = rand_exp(that.img_theta);
				auto time_to_loss = rand_exp(that.img_rho * neu.size());
				if (time_to_gain > time_to_go && time_to_loss > time_to_go) {
					break;
				}
				// gain or loss event
				/* This works because gene-gain/loss is a poisson process
				 * and as such the individual events are exponentially
				 * distributed. This means the distribution is memoryless and
				 * thus past events can be ignored.
				 */
				if (time_to_gain < time_to_loss || neu.empty()) {
					neu.emplace_back(that.loci_length, -1, locus_id++);
					acc_loci.resize(acc_loci.size() + 1);
					that.ref_acc.push_back(*(neu.end() - 1));
					time += time_to_gain;
				} else {
					assert(neu.size() != 0);
					auto loser = rand_int(0, neu.size());
					// pick one
					std::swap(neu[loser], *(neu.end() - 1));
					neu.pop_back();
					time += time_to_loss;
				}
			}

			stack.push_back(neu);
		},
		[&stack, &acc_loci,
		 img_core_size = this->img_core_size ](const tree_node &self) {
			if (self.is_branch()) {
				return;
			}

			auto &top = *(stack.end() - 1);
			auto genome_id = self.get_index();
			for (auto &loc : top) {
				loc.set_genome_id(genome_id);
				acc_loci[loc.get_locus_id() - img_core_size].push_back(loc);
			}
		},
		[&stack](const tree_node &) {
			stack.pop_back(); //
		});

	assert(stack.empty());

	// filter out empty loci which were fully lost
	assert(ref_acc.size() == acc_loci.size());
	auto it =
		remove_if_i(ref_acc.begin(), ref_acc.end(),
					[&acc_loci](size_t i) { return acc_loci.at(i).empty(); });
	ref_acc.erase(it, ref_acc.end());

	auto is_empty = [](const std::vector<locus> &set) { return set.empty(); };

	// erase-remove-idiom
	auto empty_it = std::remove_if(acc_loci.begin(), acc_loci.end(), is_empty);
	acc_loci.erase(empty_it, acc_loci.end());

	this->acc_loci = acc_loci;
}

std::vector<locus> img_model::get_reference()
{
	auto ret = ref_core;
	ret.reserve(ref_core.size() + ref_acc.size());

	std::copy(ref_acc.begin(), ref_acc.end(), std::back_inserter(ret));

	return ret;
}

std::vector<locus> img_model::get_core()
{
	return ref_core; // implicit copy
}

std::vector<locus> img_model::get_accessory()
{
	return ref_acc; // implicit copy
}

std::vector<locus> img_model::get_genome(ssize_t genome_id)
{
	auto ret = std::vector<locus>();
	auto inserter = std::back_inserter(ret);

	for (auto &locus_set : cor_loci) {
		std::copy_if(locus_set.begin(), locus_set.end(), inserter,
					 [&genome_id](const auto &loc) {
						 return loc.get_genome_id() == genome_id;
					 });
	}

	for (auto &locus_set : acc_loci) {
		std::copy_if(locus_set.begin(), locus_set.end(), inserter,
					 [&genome_id](const auto &loc) {
						 return loc.get_genome_id() == genome_id;
					 });
	}

	return ret;
}

std::vector<locus> img_model::get_locus(ssize_t some_number)
{
	if (some_number < cor_loci.size()) {
		return cor_loci.at(some_number); // throws on some_number < 0
	} else {
		return acc_loci.at(some_number - cor_loci.size()); // may throw
	}
}

std::string img_model::get_coalescent()
{
	return (coalescent.end() - 1)->to_newick();
}
