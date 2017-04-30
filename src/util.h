#pragma once
#include <algorithm>
#include <string>
#include <vector>

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

template <class ForwardIt, class UnaryPredicate>
ForwardIt find_if_i(ForwardIt first, ForwardIt last, UnaryPredicate p)
{
	auto very_first = first;
	for (; first != last; ++first) {
		if (p(std::distance(very_first, first))) {
			return first;
		}
	}
	return first;
}

template <class ForwardIt, class UnaryPredicate>
ForwardIt remove_if_i(ForwardIt first, ForwardIt last, UnaryPredicate p)
{
	auto very_first = first;
	first = find_if_i(first, last, p);
	if (first != last)
		for (ForwardIt i = first; ++i != last;)
			if (!p(std::distance(very_first, i))) *first++ = std::move(*i);
	return first;
}

// technically top could be `noexcept`.
template <typename T> T &top(std::vector<T> &v)
{
	return *(v.end() - 1);
}

template <typename T> const T &top(const std::vector<T> &v)
{
	return *(v.end() - 1);
}
