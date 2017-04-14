#include "img.h"
#include "config.h"
#include "global.h"
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

void img_model::simulate()
{
	// generate coalescent
	auto n = num_genomes;

	// Use a pool for allocating the tree nodes. Nodes 0 to n are leafes, the
	// rest are internal nodes, with the last being the root.
	auto pool = std::vector<tree_node>(2 * n - 1);
	auto &root = *(pool.end() - 1);

	for (size_t i = 0; i < n; i++) {
		pool[i] = tree_node::leaf(i);
	}

	auto indirect = std::vector<tree_node *>();
	indirect.reserve(2 * n - 1);
	for (auto &tn : pool) {
		indirect.push_back(&tn);
	}

	auto rand_int = [](size_t lower, size_t upper) {
		auto rng_help = std::uniform_int_distribution<size_t>(lower, upper - 1);
		return rng_help(RNG);
	};

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
	auto rand_exp = [](double arg) {
		auto rng_dist = std::exponential_distribution<double>(arg);
		return rng_dist(RNG);
	};

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

	std::cout << indirect[0]->to_newick() << std::endl;

	// create core sequences
	// generate pan genome and create sequences
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
