#ifndef UTILS_HPP_INCLUDED
#define UTILS_HPP_INCLUDED

#include <ctime>
#include <random>
#include <set>
#include <vector>
#include <bitset>

//=======================================================

int power(int a, int b);

//=======================================================

class Timewatch
{
public:
	clock_t timestamp;
	Timewatch();
	float watch();
};

//=======================================================

class Random
{
public:
	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution;
	Random(int i, int s);
	int next();
};

//=======================================================

Timewatch::Timewatch()
{
	timestamp = clock();
}

float Timewatch::watch()
{
	float elapsed_secs = float(clock() - timestamp) / CLOCKS_PER_SEC;
	timestamp = clock();
	return elapsed_secs;
}

//=======================================================

Random::Random(int i, int s) :
distribution(i,s)
{
}

int Random::next()
{
	return distribution(generator);
}

//=======================================================

template <size_t N>
using mm_bitset = std::bitset<N>;

//=======================================================

class mm_vector_with_properties_options
{
public:
	std::vector<std::set<int>> rnd_sets;
};

template <size_t N>
class mm_vector_with_properties
{
public:
	mm_bitset<N> v;
	int bit_count;
	static const int rnd_size = 0;
	static const int rnd_count = 8;
	mm_bitset<rnd_size> r;
	void calculate_properties(mm_vector_with_properties_options & o);
	static void make_options(mm_vector_with_properties_options & o);
};

template <size_t N>
void mm_vector_with_properties<N>::calculate_properties(mm_vector_with_properties_options & o)
{
	bit_count = v.count();
	for (int i = 0; i < rnd_size; ++i)
	{
		r.set(i,false);
		for (std::set<int>::iterator j = o.rnd_sets[i].begin(); j != o.rnd_sets[i].end(); ++j)
		{
			r[i] = r[i] != v[*j];
		}
	}
}

template <size_t N>
void mm_vector_with_properties<N>::make_options(mm_vector_with_properties_options & o)
{
	Random g(0, N-1);
	for (int i = 0; i < rnd_size; ++i)
	{
		std::set<int> c;
		for (int j = 0; j < rnd_count; ++j)
			{
				int n = g.next();
				while (c.find(n) != c.end())
					n = g.next();
				c.insert(n);
			}
		o.rnd_sets.push_back(c);
	}
}

//=======================================================

template <size_t N>
class Vectors_Presolve_Data
{
public:
	mm_bitset<N> mor;
	Vectors_Presolve_Data();
	void add_vector(const mm_bitset<N> & a);
	bool check(const mm_bitset<N> & a);
};

template <size_t N>
Vectors_Presolve_Data<N>::Vectors_Presolve_Data()
{
	mor.reset();
}

template <size_t N>
void Vectors_Presolve_Data<N>::add_vector(const mm_bitset<N> & a)
{
	mor |= a;
}

template <size_t N>
bool Vectors_Presolve_Data<N>::check(const mm_bitset<N> & a)
{
	if ((~mor & a).any())
		return false;
	return true;
}

//=======================================================

int power(int a, int b)
{
    int result = a;
    while (--b)
        result *= a;
    return result;
}

//=======================================================

#endif // UTILS_HPP_INCLUDED
