#include "img.h"
#include "config.h"
#include "global.h"
#include <iostream>
#include <sstream>
#include <string>

class tree_node
{
	tree_node *parent = nullptr;
	tree_node *left_child = nullptr;
	tree_node *right_child = nullptr;
	double up_time = 0.0;
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

	bool is_branch() const
	{
		return left_child != nullptr;
	}

	bool is_leaf() const
	{
		return left_child == nullptr;
	}

	ssize_t get_index() const
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
		};

		traverse(pre, process, post);

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
	auto pool = std::vector<tree_node>(2 * n - 1);

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

	std::cerr << "n: " << n << std::endl;
	for (size_t i = 0; i < n - 1; i++) {
		auto p = n + i;
		std::cerr << "i: " << i << " p: " << p << std::endl;
		auto c = rand_int(0, n - i);
		std::cerr << "c: " << c << std::endl;
		indirect.at(p)->set_left(indirect[c]); // also sets [c].parent
		indirect[c] = indirect[n - i - 1];

		auto d = rand_int(0, n - i - 1);
		std::cerr << "d: " << d << std::endl;
		indirect[p]->set_right(indirect[d]);
		indirect[d] = indirect[p];
	}

	std::cerr << indirect[0]->to_newick() << std::endl;

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
