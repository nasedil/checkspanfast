/**
 * @file utils.hpp:  various utilities for
 * Cube Product Checker program.
 *
 * Contains power function,
 * Timewatch class,
 * Random class,
 * mm_bitset type alias,
 * mm_vector_with_properties class,
 * mm_vector_with_properties_options class,
 * Vectors_Presolve_Data class.
 *
 * @author Eugene Petkevich
 * @version pre-alpha
 */

#ifndef UTILS_HPP_INCLUDED
#define UTILS_HPP_INCLUDED

#include <ctime>
#include <random>
#include <set>
#include <vector>
#include <bitset>
#include <chrono>

using namespace std;

//=============================================================================

/**
 * Type alias for bitset.
 */
template <size_t N>
using mm_bitset = bitset<N>;

//=============================================================================

int power(int a, int b); /// Rise a to power of b.
int popcount(int x); /// Return number of set bits in a number.

//=============================================================================

/**
 * Timewatch class, to measure time intervals.
 */
class Timewatch {
public:
    double timestamp; /// Timestamp of the last call.
    Timewatch();
    double watch(); /// Get interval duration since last call.
};

//=============================================================================

/**
 * Random class, for integer uniform distribution.
 */
class Random {
public:
    mt19937 generator;
    uniform_int_distribution<int>* distribution;
    Random();
    ~Random();
    void init(int i, int s); /// Init with interval [i,s].
    int next(); /// Return random integer form the interval.
};

//=============================================================================

/**
 * Options for adding randomized bits to vectors.
 *
 * Contains sets of bit indices to add up
 * for making additional bits in vectors.
 */
class mm_vector_with_properties_options
{
public:
    vector<set<int>> rnd_sets; /// Sets of the bit indices.
};

/**
 * Bitset vector with additional properties.
 *
 * Additional bits are added
 * which are linear combinations of original bits.
 *
 * @param N: size of the bitsets.
 */
template <size_t N>
class mm_vector_with_properties {
public:
    mm_bitset<N> v; /// Main bitset.
    int v_count; /// number of set bits in v
    static const int rnd_size = 0; /// Number of terms for linear combination.
    static const int rnd_count = 0; /// Number of additional bits.
    mm_bitset<rnd_size> r; /// Additional bits.
    int r_count; /// number of set bits in r
    void calculate_properties(mm_vector_with_properties_options& o); /// Calculate additional bits.
    static void make_options(mm_vector_with_properties_options& o); /// Make options randomly.
};

//=============================================================================

/**
 * Bitset set properties.
 *
 * Stores some properties of a set of bitsets,
 * to do fast check for linear independence.
 *
 * @param N: size of the bitsets.
 */
template <size_t N>
class Vectors_Presolve_Data
{
public:
    mm_bitset<N> mor; /// Arithmetical OR of the all bitsets.
    Vectors_Presolve_Data();
    void add_vector(const mm_bitset<N>& a); /// Add a bitset.
    void presolve(); /// presolve
    bool check(const mm_bitset<N>& a) const; /// Check if a vector could be in the span.
};

//=============================================================================
//=============================================================================

/**
 * Constructor.
 */
Timewatch::Timewatch()
{
    timestamp = omp_get_wtime();
}

/**
 * Get interval duration since last call.
 *
 * @return interval duration in seconds.
 */
double
Timewatch::watch()
{
    double elapsed_secs = double(omp_get_wtime() - timestamp);
    timestamp = omp_get_wtime();
    return elapsed_secs;
}

//=============================================================================

Random::Random() :
    generator(chrono::system_clock::now().time_since_epoch().count())
{
    distribution = new uniform_int_distribution<int>(0, 1);
}

Random::~Random()
{
    delete distribution;
}

/**
 * Init with interval [i,s].
 *
 * @param lower bound of the interval (inclusive);
 * @param upper bound of the interval (inclusive).
 */
void
Random::
init(int i, int s)
{
    delete distribution;
    distribution = new uniform_int_distribution<int>(i, s);
}

/**
 * Return random integer form the interval.
 *
 * @return resulting random number.
 */
int
Random::
next()
{
    return (*distribution)(generator);
}

//=============================================================================

/**
 * Calculate properties of the vector.
 *
 * @param o: options.
 */
template <size_t N>
void
mm_vector_with_properties<N>::
calculate_properties(mm_vector_with_properties_options& o)
{
    v_count = v.count();
    for (int i = 0; i < rnd_size; ++i) {
        r.set(i,false);
        for (set<int>::iterator j = o.rnd_sets[i].begin(); j != o.rnd_sets[i].end(); ++j) {
            r[i] = r[i] != v[*j];
        }
    }
    r_count = r.count();
}

/**
 * Generate options randomly.
 *
 * @param o: options that will be changed.
 */
template <size_t N>
void
mm_vector_with_properties<N>::
make_options(mm_vector_with_properties_options& o)
{
    Random g;
    g.init(0, N-1);
    for (int i = 0; i < rnd_size; ++i) {
        set<int> c;
        for (int j = 0; j < rnd_count; ++j) {
            int n = g.next();
            while (c.find(n) != c.end()) {
                n = g.next();
            }
            c.insert(n);
        }
        o.rnd_sets.push_back(c);
    }
}

//=============================================================================

/**
 * Constructor.
 */
template <size_t N>
Vectors_Presolve_Data<N>::
Vectors_Presolve_Data()
{
    mor.reset();
}

/**
 * Add a bitvector to the current set.
 *
 * @param a: bitvector to add.
 */
template <size_t N>
void
Vectors_Presolve_Data<N>::
add_vector(const mm_bitset<N>& a)
{
    mor |= a;
}

/**
 * Presolve for the current set.
 */
template <size_t N>
void
Vectors_Presolve_Data<N>::
presolve()
{
    mor = ~mor;
}

/**
 * Check if a vector could be in the span of current set.
 *
 * @param a: a vector to check.
 *
 * @return false if the vector is definitely not in the span.
 */
template <size_t N>
bool Vectors_Presolve_Data<N>::check(const mm_bitset<N>& a) const
{
    if ((mor & a).any())
        return false;
    return true;
}

//=============================================================================

/**
 * Rise a to power of b.
 *
 * @param a: base;
 * @param b: exponent.
 *
 * @return a^b.
 */
int
power(int a, int b)
{
    if (b == 0) {
        return 1;
    }
    int result = a;
    while (--b) {
        result *= a;
    }
    return result;
}

//=============================================================================

/**
 * Return number of set bits in a number (max 16 bit).
 *
 * @param x: the number;
 *
 * @return the number of set bits.
 */
int popcount(int x)
{
    bitset<16> b(x);
    return b.count();
}

//=============================================================================

#endif // UTILS_HPP_INCLUDED
