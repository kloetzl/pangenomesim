#pragma once
#include <string>
#include <vector>

class gene
{
	std::string nucl = {};
	std::string m_name = {};
	ssize_t genome_id = -1;
	ssize_t gene_id = -1;

  public:
	gene() = default;
	gene(std::string, ssize_t, ssize_t);
	gene(std::string, std::string, ssize_t, ssize_t);
	gene(ssize_t, ssize_t, ssize_t);

	std::string name() const;
	const std::string &get_nucl() const;
	ssize_t get_genome_id() const;
	ssize_t get_gene_id() const;
	void set_genome_id(ssize_t);

	gene mutate(double) const;
	static gene generate_random();
	static std::vector<gene> vector_mutate(const std::vector<gene> &, double);

	std::string to_fasta() const
	{
		auto ret = std::string(">");
		ret += name();
		ret += "\n";
		ret += nucl;
		ret += "\n";
		return ret;
	}
};
