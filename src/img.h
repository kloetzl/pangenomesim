#pragma once
#include "locus.h"
#include "util.h"
#include <cassert>
#include <string>
#include <vector>

class tree_node;

class img_model
{
  protected:
	size_t loci_length = 100;
	size_t num_genomes = 3;
	size_t seed = 1729;
	double img_theta = 0.1;
	double img_rho = 0.1;
	size_t img_core_size = 4;
	double mut_rate = 0.01;
	std::vector<locus> ref_core = {};
	std::vector<locus> ref_acc = {};
	std::vector<std::vector<locus>> cor_loci = {};
	std::vector<std::vector<locus>> acc_loci = {};
	std::vector<tree_node> coalescent = {};
	std::vector<std::vector<double>> distmatrix = {};

  public:
	img_model() = default;

	void parse_param(std::string, std::string);
	std::string parameters() const;
	void simulate();

	std::vector<locus> get_reference();
	std::vector<locus> get_core();
	std::vector<locus> get_accessory();
	std::vector<locus> get_genome(ssize_t);
	std::vector<locus> get_locus(ssize_t);

	std::string get_coalescent();
	const std::vector<std::vector<double>> &get_distmatrix() const noexcept
	{
		return distmatrix;
	}

	size_t get_num_loci() const noexcept
	{
		return cor_loci.size() + acc_loci.size();
	}

	size_t get_num_genomes() const noexcept
	{
		return num_genomes;
	}

	auto get_loci_length() const noexcept
	{
		return loci_length;
	}

  private:
	std::vector<locus> seq_from_root(const tree_node &root, size_t locus_id);
};

class tree_node
{
	//< @brief a pointer to the parent, if existing.
	tree_node *parent = nullptr;
	tree_node *left_child = nullptr;
	tree_node *right_child = nullptr;
	//< @brief The length of the branch from the parent.
	double up_time = 0.0;
	//< The absolute time
	double abs_time = 0.0;
	ssize_t index = -1;

  public:
	/**
	 * @brief Create a leaf.
	 * @param idx - The index of the leaf.
	 * @returns a new leaf
	 */
	static tree_node leaf(ssize_t idx)
	{
		auto ret = tree_node();
		ret.index = idx;
		return ret;
	}
	tree_node() = default;
	~tree_node() = default;

	/**
	 * @brief Set the left child.
	 * @param new_left - The new left child.
	 */
	void set_left(tree_node *new_left)
	{
		left_child = new_left;
		if (left_child) {
			left_child->parent = this;
		}
	}

	/**
	 * @brief Set the right child.
	 * @param new_right - The new right child.
	 */
	void set_right(tree_node *new_right)
	{
		right_child = new_right;
		if (right_child) {
			right_child->parent = this;
		}
	}

	/**
	 * @brief Set the time relative to parent.
	 * @param t - The new relative time.
	 */
	void set_time(double t)
	{
		assert(t > 0);
		up_time = t;
	}

	/**
	 * @brief Set the absolute time. (Now is 0.)
	 * @param t - The new absolute time.
	 */
	void set_abs_time(double t)
	{
		assert(t > 0);
		abs_time = t;
	}

	/**
	 * Compute the branch length.
	 */
	void compute_rel_time()
	{
		assert(has_parent());

		auto branch_length = get_parent().get_abs_time() - get_abs_time();
		set_time(branch_length);
	}

	/**
	 * @brief Get the realtive time.
	 * @returns the time relative to parent.
	 */
	double get_time() const noexcept
	{
		return up_time;
	}

	/**
	 * @brief get the absolute time
	 * @returns the absolute time.
	 */
	double get_abs_time() const noexcept
	{
		return abs_time;
	}

	/**
	 * @brief Check whether the current node  has a parent node. This is only
	 * false for the root.
	 * @returns true iff the current node has a parent.
	 */
	bool has_parent() const noexcept
	{
		return parent != nullptr;
	}

	/**
	 * @brief Get a reference to the parent. Must only be called on node that
	 * have a parent.
	 * @returns the parent node.
	 */
	tree_node &get_parent() noexcept
	{
		return *parent;
	}

	/**
	 * @brief Get a reference to the parent. Must only be called on node that
	 * have a parent.
	 * @returns the parent node.
	 */
	const tree_node &get_parent() const noexcept
	{
		return *parent;
	}

	/**
	 * @brief Checks whether this is an internal node.
	 * @returns true iff the current node is branching.
	 */
	bool is_branch() const noexcept
	{
		return left_child != nullptr;
	}

	/**
	 * @brief Checks whether the current node is a leaf.
	 * @returns true iff the current node is a leaf.
	 */
	bool is_leaf() const noexcept
	{
		return left_child == nullptr;
	}

	/**
	 * @brief Get the index of the node. For leaves, this is a positive number.
	 * For branches it is -1.
	 * @returns the index.
	 */
	ssize_t get_index() const noexcept
	{
		return index;
	}

	/**
	 * @brief Execute an in-order traversal of the tree. Call the passed
	 * function on each node.
	 *
	 * @param process - A function to execute on every node.
	 */
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

	/**
	 * @brief Execute an in-order traversal of the tree. Call the passed
	 * function on each node.
	 *
	 * @param process - A function to execute on every node.
	 */
	template <typename Func> void traverse(const Func &process) const
	{
		if (left_child) {
			left_child->traverse(process);
		}
		process(*this);
		if (right_child) {
			right_child->traverse(process);
		}
	}

	/**
	 * @brief Execute an in-order traversal of the tree. Three function may be
	 * passed executed on entry, processing, and exit of a node, respectively.
	 *
	 * @param pre - A function to be executed whenever the traversal enters a
	 * node.
	 * @param process - A function to be executed on the node. In case of a
	 * branch it is executed in between the left and right branch.
	 * @param post - A function to be executed when the traversal leaves a node.
	 */
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

	/**
	 * @brief Execute an in-order traversal of the tree. Three function may be
	 * passed executed on entry, processing, and exit of a node, respectively.
	 *
	 * @param pre - A function to be executed whenever the traversal enters a
	 * node.
	 * @param process - A function to be executed on the node. In case of a
	 * branch it is executed in between the left and right branch.
	 * @param post - A function to be executed when the traversal leaves a node.
	 */
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

	/**
	 * @brief Convert tree to newick format.
	 *
	 * @returns String in Newick format.
	 */
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
				ret += genome_name(self.get_index());
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
