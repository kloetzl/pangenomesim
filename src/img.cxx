#include "img.h"
#include "config.h"
#include "global.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>

class tree_node
{
	tree_node *parent = nullptr;
	tree_node *left_child = nullptr;
	tree_node *right_child = nullptr;
	double up_time = 0.0;
	double abs_time = 0.0;
	ssize_t index = -1; //?

  public:
	static tree_node leaf(ssize_t idx)
	{
		auto ret = tree_node();
		ret.index = idx;
		return ret;
	}
	tree_node() = default;
	~tree_node() = default;

	void set_left(tree_node *new_left)
	{
		left_child = new_left;
		left_child->parent = this;
	}

	void set_right(tree_node *new_right)
	{
		right_child = new_right;
		right_child->parent = this;
	}

	void set_time(double t)
	{
		assert(t > 0);
		up_time = t;
	}

	void set_abs_time(double t)
	{
		assert(t > 0);
		abs_time = t;
	}

	double get_time() const noexcept
	{
		return up_time;
	}

	double get_abs_time() const noexcept
	{
		return abs_time;
	}

	bool has_parent() const noexcept
	{
		return parent != nullptr;
	}

	tree_node &get_parent() noexcept
	{
		return *parent;
	}

	bool is_branch() const noexcept
	{
		return left_child != nullptr;
	}

	bool is_leaf() const noexcept
	{
		return left_child == nullptr;
	}

	ssize_t get_index() const noexcept
	{
		return index;
	}

	template <typename Func> void traverse(const Func &process)
	{
		if (left_child) {
			left_child->traverse(process);
		}
		process(*this);
		if (right_child) {
			right_child->traverse(process);
		}
	}

	template <typename Func1, typename Func2, typename Func3>
	void traverse(const Func1 &pre, const Func2 &process, const Func3 &post)
	{
		pre(*this);
		if (left_child) {
			left_child->traverse(pre, process, post);
		}
		process(*this);
		if (right_child) {
			right_child->traverse(pre, process, post);
		}
		post(*this);
	}

	template <typename Func1, typename Func2, typename Func3>
	void traverse(const Func1 &pre, const Func2 &process,
				  const Func3 &post) const
	{
		pre(*this);
		if (left_child) {
			left_child->traverse(pre, process, post);
		}
		process(*this);
		if (right_child) {
			right_child->traverse(pre, process, post);
		}
		post(*this);
	}

	std::string to_newick()
	{
		auto ret = std::string();

		auto pre = [&ret](const tree_node &self) {
			if (self.is_branch()) {
				ret += "(";
			}
		};
		auto process = [&ret](const tree_node &self) {
			if (self.is_leaf()) {
				ret += std::to_string(self.get_index());
			} else {
				ret += ",";
			}
		};
		auto post = [&ret](const tree_node &self) {
			if (self.is_branch()) {
				ret += ")";
			}
			if (self.up_time != 0.0) {
				ret += ":";
				ret += std::to_string(self.up_time);
			}
		};

		traverse(pre, process, post);

		ret += ";";
		return ret;
	}
};

void img_model::parse_param(std::string key, std::string value)
{
	if (key == "img_theta") {
		this->img_theta = std::stoul(value);
		return;
	}
	if (key == "img_rho") {
		this->img_rho = std::stoul(value);
		return;
	}
	evo_model::parse_param(key, value);
}

std::string img_model::parameters() const
{
	auto str = std::stringstream();

	str << "pangenomesim " << VERSION << "\n\n";
	str << "full options:\n";
	str << "--param gene_length=" << loci_length << "\n";
	str << "--param num_loci=" << num_loci << "\n";
	str << "--param num_genomes=" << num_genomes << "\n";
	str << "--param img_theta=" << img_theta << "\n";
	str << "--param img_rho=" << img_rho << "\n";
	str << "--param seed=" << seed << "\n";

	return str.str();
}

std::vector<tree_node> coalescent(size_t n)
{
	// Use a pool for allocating the tree nodes. Nodes 0 to n are leaves, the
	// rest are internal nodes, with the last being the root.
	auto pool = std::vector<tree_node>(2 * n - 1);
	auto &root = *(pool.end() - 1);

	generate_i(pool.begin(), n, [](size_t i) {
		return tree_node::leaf(i); //
	});

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

	// compute relative times
	root.traverse([](tree_node &self) {
		if (self.has_parent()) {
			auto branch_length =
				self.get_parent().get_abs_time() - self.get_abs_time();
			self.set_time(branch_length);
		}
	});

	return pool;
}

std::vector<locus> seq_from_root(const tree_node &root, size_t sample_size,
								 size_t loci_length, double rate,
								 size_t locus_id)
{
	// create core sequences
	auto leaves = std::vector<locus>(sample_size);
	auto stack = std::vector<locus>();
	stack.reserve(sample_size);
	auto seq = locus(loci_length, -1, locus_id); // random root seq
	stack.push_back(seq);
	// ref_core.push_back(seq);
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

auto locus_set_mutate(const std::vector<locus> &set, double rate)
{
	auto ret = std::vector<locus>();
	ret.reserve(set.size());
	for (const auto &loc: set){
		ret.push_back(loc.mutate(rate));
	}
	return ret;
}

void img_model::simulate()
{
	// generate coalescent
	auto pool = coalescent(num_genomes);
	auto &root = *(pool.end() - 1);

	std::cout << root.to_newick() << std::endl;

	auto rate = 0.1;

	// create core sequences
	generate_i(std::back_inserter(loci), img_core_size, [&](size_t locus_id) {
		return seq_from_root(root, num_genomes, loci_length, rate, locus_id);
	});

	using loci_set = std::vector<locus>;

	// generate pan genome and create sequences
	auto acc_loci = std::vector<loci_set>(num_genomes);
	auto acc_locus = std::vector<locus>();

	auto start = std::vector<locus>();
	auto start_size = rand_poisson(img_theta / img_rho);
	std::cerr << "acc start size: " << start_size << std::endl;
	start.reserve(start_size);
	generate_i(std::back_inserter(start), start_size, [=](size_t locus_id) {
		return locus(loci_length, -1, locus_id);
	});

	auto locus_counter = start_size;
	ref_acc = start; // copy

	auto stack = std::vector<loci_set>();
	stack.reserve(num_genomes);
	stack.push_back(start);

	root.traverse([&stack, &rate, &that = *this,
				   &locus_id = locus_counter ](const tree_node &self) {
		if (!self.has_parent()) {
			return; // root
		}

		const auto &top = *(stack.end() - 1);
		// simulate evolution
		auto neu = locus_set_mutate(top, rate * self.get_time());
		auto time = 0.0;
		while (time < self.get_time()) {
			std::cerr << time << " " << self.get_time() << std::endl;
			auto time_to_go = self.get_time() - time;
			auto time_to_gain = exp(that.img_theta);
			auto time_to_loss = exp(that.img_rho * neu.size());
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
				std::cerr << "gain!" << std::endl;
				neu.emplace_back(that.loci_length, -1, locus_id++);
				that.ref_acc.push_back(*(neu.end() - 1));
				time += time_to_gain;
			} else {
				std::cerr << "loss!" << std::endl;
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
				  [&stack, &acc_loci](const tree_node &self) {
					  if (self.is_branch()) {
						  return;
					  }

					  auto &top = *(stack.end() - 1);
					  auto index = self.get_index();
					  for (auto &loc : top) {
						  loc.set_genome_id(index);
					  }

					  acc_loci[index] = top;
				  },
				  [&stack](const tree_node &) {
					  stack.pop_back(); //
				  });

	assert(stack.empty());
	acc = acc_loci;
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
