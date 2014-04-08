#ifndef UTILS_HPP_INCLUDED
#define UTILS_HPP_INCLUDED

#include <ctime>
#include <random>

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

template <size_t N>
class mm_vector_with_properties
	{
	public:
		mm_bitset<N> v;
		int bit_count;
		const int mask_count = N/32;
		//int mask0[NM/32];
		//int mask1[NM/32];

		void calculate_properties()
		{
			bit_count = v.count();
			/*
			for (int i = 0; i < mask_count; ++i)
			{

			}
			*/
		}
	};

//=======================================================

#endif // UTILS_HPP_INCLUDED
