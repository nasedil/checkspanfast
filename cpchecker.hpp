/**
 * @file mmchecker.hpp:  Cube_Product_Checker template class
 * and its implementation.
 *
 * Class is for finding Strassen-like algorithms for cube product.
 *
 * @author Eugene Petkevich
 * @version pre-alpha
 */

#ifndef MMCHECKER_HPP_INCLUDED
#define MMCHECKER_HPP_INCLUDED

#include <algorithm>
#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <map>
#include <string>
#include <numeric>
#include <omp.h>

#include <boost/dynamic_bitset.hpp>
#include <boost/algorithm/string/join.hpp>

#include "slae.hpp"
#include "utils.hpp"

using namespace std;

//=============================================================================

class Solution_Properties;

//=============================================================================
//=============================================================================

/**
 * Main class, calculates all nesessary information
 * for given template parameters:
 *
 * @param N: size of the cube;
 * @param D: dimension of the cube;
 * @param NM: number of bits in multiplication vectors ((N^D)^D);
 * @param NMH: number of elements in the array (N^D).
 */
template <int N, int D, size_t NM, size_t NMH>
class Cube_Product_Checker {
public:
    typedef mm_bitset<NM> Multiplication_Vector;
    typedef mm_bitset<NMH> Multiplication_Part_Vector;
    typedef set<int> Candidate;
    int length; /// size of cube (N)
    int dimension; /// dimension of cube (D)
    int element_count; /// number of elements in cube (N^D)
    mm_vector_with_properties_options vector_options; /// options for vector properties
    mm_vector_with_properties<NM>* r_vectors; /// result vectors
    mm_vector_with_properties<NM>* m_vectors; /// multiplication vectors
    bool owns_arrays; /// If arrays are created by the object.
    int m_count; /// number of non-zero multiplication vectors = (2^(N^D)-1)^D
    int m_length; /// number non-zero element sums = 2^(N^D)-1
    int f_count; /// number of vectors to choose for a basis (for minimal improvement it is N^(D+1)-1)
    set<int> good_vectors_indexes; /// set of indexes for vectors that are in the current span
    set<int> n_vectors_indexes; /// set of indexes of current chosen vectors
    set<Candidate> neighbours; /// set of neighbours;
    set<set<int>> solutions; /// set of found unique solutions
    map<set<int>, int> solution_distribution; /// set of all found solutions
    int iteration_count; /// number of finished iterations
    int local_max_iterations; /// number of iterations in probable local maximum
    int restarts; /// number of restarts
    int local_iterations; /// number of iterations
    int checked_sets_count; /// number of sets checked
    int raw_solution_count; /// number of found solutions, including duplicates
    int best_result; /// maximum number of vectors in the span found
    int lin_dependent_sets;
    Timewatch tw; /// timer that is used for getting calculation time
    Random rnd; /// random number generator
    int thread_number; /// thread number
    map<Candidate, int> candidate_cache; /// cache of already checked candidates
    int cache_limit; /// limit of the cache size
    int cache_hits; /// cache hits
    bool* stop_signal; /// if execution should be stopped
    int bit_check_hits; /// number of bit check hits
    int gaussian_eliminations; /// number of gaussian eliminations
    set<set<int>> top_best_solutions; /// best solutions found so far
    int top_best_result; /// best result so far
    set<set<int>> candidate_space; /// candidate space

    //=============--- constructors and destructors
    Cube_Product_Checker();
    ~Cube_Product_Checker();

    //=============--- index operations
    /// return linear index in a bit vector by its indexes in matrices
    int get_vector_index(int ai, int aj, int bi, int bj) const; /// for 2-dimensional case
    int get_vector_index(int ai, int aj, int ak, int bi, int bj, int bk, int ci, int cj, int ck) const; /// for 3-d case
    /// return linear index in a bit vector by combined indexes in matrices
    int get_vector_index(int a, int b) const; /// for 2-d case
    int get_vector_index(int a, int b, int c) const; /// for 3-d case
    /// return indices from bit index
    void decode_indices_from_index(int index, int& ai, int& aj, int& bi, int& bj) const; /// for 2-d case
    void decode_indices_from_index(int index, int& ai, int& aj, int& ak, int& bi, int& bj, int& bk, int& ci, int& cj, int& ck) const; /// for 3-d case
    /// return linear index of an element in a matrix
    int get_element_index(int i, int j) const; /// for 2-d case
    int get_element_index(int i, int j, int k) const; /// for 3-d case
    /// return index of multiplication vector in the set
    int get_m_index(int i, int j) const; /// for 2-d case
    int get_m_index(int i, int j, int k) const; /// for 3-d case
    vector<int> decode_m_index(int index) const; /// return the coefficients from m-vector index

    //=============--- Initial calculations
    void init(int mult_count, const string& filename, const string& space_filename); /// calculate all properties
    void init(const Cube_Product_Checker& cpc); /// Link to all properties in other object.
    void calculate_r_vectors(); /// write result vectors to array
    void calculate_m_vectors(); /// write multiplication vectors to array

    //=============--- checking routines
    void clear_sets(); /// clear current sets of vectors
    void add_vector_to_set(int index); /// add vector to current sets
    bool check_vectors_for_goodness(); /// check current set of vectors
    void clear_statistics(); /// clear statistics of solutions
    void make_random_candidate(); /// make random candidate solution
    void make_candidate_from_space(); /// make random candidate from restricted space
    bool check_cache(); /// check if candidate result is in cache
    void update_cache(); /// update the cache with new result

    //=============--- searching for solution
    bool check_for_good_vectors(); /// check all solution space
    bool check_for_good_vectors_randomized(); /// do random search
    bool solve_hill_climbing(int local_max_limit, bool use_space); /// do local search
    void start_neighbourhood(); /// make list of neighbours

    //=============--- utilities
    void output_vector(Multiplication_Vector v) const; /// output vector to screen
    void output_vector_text(Multiplication_Vector v) const; /// output vector to screen in letters
    void save_random_samples(int size, const char* filename) const; /// save random sets to a file (for testing later)
    void read_samples_and_check(const char* filename, const char* filenameout) const; /// check sets from a file
    void output_current_state() const; // output current state of the checker
    void read_m_vectors(const string& filename); /// read m_vectors from file
    void read_candidate_space(const string& filename); /// read candidate space from file
    void write_m_vectors(const string& filename); /// write m_vectors into file

    //=============--- Statistics and results
    void save_results(const char* filename); /// save results to a file
    bool check_solution(set<int> s, Solution_Properties& sp); /// check if a solution is valid
    void save_solution_properties(const Solution_Properties& sp, const char* filename); /// output solution properties to a file
    vector<int> sum_operations_cube(int index); /// number of summation operations inside cubes
};

//=============================================================================

class Solution_Properties {
public:
    set<int> multiplication_vectors; /// multiplication vector indices
    vector<boost::dynamic_bitset<>> coefficients; /// result vector coefficients
    int operation_count; /// number of addition operations used for calculating result matrix overall
};

//=============================================================================
//=============================================================================

/**
 * Class constructor.
 *
 * Create an object for work with D-dimensional cubes of size N.
 */
template <int N, int D, size_t NM, size_t NMH>
Cube_Product_Checker<N, D, NM, NMH>::
Cube_Product_Checker() :
    owns_arrays(false)
{
    length = N;
    dimension = D;
    element_count = power(length, dimension);
    f_count = power(length,dimension+1)-1;
    m_length = power(2,element_count)-1;
    m_count = power(m_length,dimension);
    cache_limit = 3000000;
    cache_hits = 0;
    bit_check_hits = 0;
    gaussian_eliminations = 0;
    rnd.init(0, m_count-1);
#ifdef VERBOSE_OUTPUT
    cout << "Cube Product Checker has been created." << endl;
#endif // VERBOSE_OUTPUT
}

//=============================================================================

/**
 * Destructor.
 */
template <int N, int D, size_t NM, size_t NMH>
Cube_Product_Checker<N, D, NM, NMH>::
~Cube_Product_Checker()
{
    if (owns_arrays) {
        delete [] r_vectors;
        delete [] m_vectors;
        delete stop_signal;
    }
}

//=============================================================================

/**
 * Initialize the object with all necessary properties.
 *
 * @param mult_count:  number of multiplications in resulting algorithm.
 */
template <int N, int D, size_t NM, size_t NMH>
void
Cube_Product_Checker<N, D, NM, NMH>::
init(int mult_count, const string& filename, const string& space_filename)
{
    f_count = mult_count;
    stop_signal = new bool(false);
    r_vectors = new mm_vector_with_properties<NM>[element_count];
    m_vectors = new mm_vector_with_properties<NM>[m_count];
    mm_vector_with_properties<NM>::make_options(vector_options);
    owns_arrays = true;
#ifdef VERBOSE_OUTPUT
    tw.watch();
#endif // VERBOSE_OUTPUT
    calculate_r_vectors();
    for (int i = 0; i < element_count; ++i) {
        r_vectors[i].calculate_properties(vector_options);
    }
#ifdef VERBOSE_OUTPUT
    cout << "[" << tw.watch() << " s] Result vectors calculated." << endl;
#endif // VERBOSE_OUTPUT
    if (filename.length() > 0) {
        read_m_vectors(filename);
    } else {
        calculate_m_vectors();
        for (int i = 0; i < m_count; ++i) {
            m_vectors[i].calculate_properties(vector_options);
        }
    }
    if (space_filename.length() > 0) {
        read_candidate_space(space_filename);
    }
#ifdef VERBOSE_OUTPUT
    cout << "[" << tw.watch() << " s] Multiplication vectors calculated." << endl;
#endif // VERBOSE_OUTPUT
}

/**
 * Initialize the object with all necessary properties from other object.
 *
 * Link to all big data fields in the other object.
 *
 * @param cpc: the object to get main arrays from.
 */
template <int N, int D, size_t NM, size_t NMH>
void
Cube_Product_Checker<N, D, NM, NMH>::
init(const Cube_Product_Checker& cpc)
{
    f_count = cpc.f_count;
    r_vectors = cpc.r_vectors;
    m_vectors = cpc.m_vectors;
    vector_options = cpc.vector_options;
    stop_signal = cpc.stop_signal;
    candidate_space = cpc.candidate_space;
}

//=============================================================================

/**
 * Get bit index in a multiplication vector (2-dimensional case).
 *
 * @param ai: element's first index in the first matrix;
 * @param aj: element's second index in the first matrix;
 * @param bi: element's first index in the second matrix;
 * @param bj: element's second index in the second matrix;
 *
 * @return bit index.
 */
template <int N, int D, size_t NM, size_t NMH>
inline int
Cube_Product_Checker<N, D, NM, NMH>::
get_vector_index(int ai, int aj,
                 int bi, int bj) const
{
    int result = ai;
    result *= length;
    result += aj;
    result *= length;

    result += bi;
    result *= length;
    result += bj;
    return result;
}

/**
 * Get bit index in a multiplication vector (3-dimensional case).
 *
 * @param ai: element's first index in the first cube;
 * @param aj: element's second index in the first cube;
 * @param ak: element's third index in the first cube;
 * @param bi: element's first index in the second cube;
 * @param bj: element's second index in the second cube;
 * @param bk: element's third index in the second cube;
 * @param ci: element's first index in the third cube;
 * @param cj: element's second index in the third cube;
 * @param ck: element's third index in the third cube.
 *
 * @return bit index.
 */
template <int N, int D, size_t NM, size_t NMH>
inline int
Cube_Product_Checker<N, D, NM, NMH>::
get_vector_index(int ai, int aj, int ak,
                 int bi, int bj, int bk,
                 int ci, int cj, int ck) const
{
    int result = ai;
    result *= length;
    result += aj;
    result *= length;
    result += ak;
    result *= length;

    result += bi;
    result *= length;
    result += bj;
    result *= length;
    result += bk;
    result *= length;

    result += ci;
    result *= length;
    result += cj;
    result *= length;
    result += ck;
    return result;
}

//=============================================================================

/**
 * Get bit index in a multiplication vector (2-dimensional case).
 *
 * @param a: element'S index in the first matrix;
 * @param b: element'S index in the second matrix.
 *
 * @return bit index.
 */
template <int N, int D, size_t NM, size_t NMH>
inline int
Cube_Product_Checker<N, D, NM, NMH>::
get_vector_index(int a, int b) const
{
    int result = a;
    result *= element_count;
    result += b;
    return result;
}

/**
 * Get bit index in a multiplication vector (3-dimensional case).
 *
 * @param a: element's index in the first cube;
 * @param b: element's index in the second cube;
 * @param c: element's index in the third cube.
 *
 * @return bit index.
 */
template <int N, int D, size_t NM, size_t NMH>
inline int
Cube_Product_Checker<N, D, NM, NMH>::
get_vector_index(int a, int b, int c) const
{
    int result = a;
    result *= element_count;
    result += b;
    result *= element_count;
    result += c;
    return result;
}

//=============================================================================

/**
 * Get element indices in matrices from bit index in the multiplication vector
 * (2-dimensional case).
 *
 * @param index: bit index in a multiplication vector;
 *
 * @param ai: element's first index in the first matrix;
 * @param aj: element's second index in the first matrix;
 * @param bi: element's first index in the second matrix;
 * @param bj: element's second index in the second matrix.
 */
template <int N, int D, size_t NM, size_t NMH>
inline void
Cube_Product_Checker<N, D, NM, NMH>::
decode_indices_from_index(int index,
                          int& ai, int& aj,
                          int& bi, int& bj) const
{
    bj = index % length;
    index /= length;
    bi = index % length;;
    index /= length;
    aj = index % length;
    index /= length;
    ai = index;
    return;
}

/**
 * Get element indices in cubes from bit index in a multiplication vector
 * (2-dimensional case).
 *
 * @param index: bit index in the multiplication vector;
 *
 * @param ai: element's first index in the first cube;
 * @param aj: element's second index in the first cube;
 * @param ak: element's third index in the first cube;
 * @param bi: element's first index in the second cube;
 * @param bj: element's second index in the second cube;
 * @param bk: element's third index in the second cube;
 * @param ci: element's first index in the third cube;
 * @param cj: element's second index in the third cube;
 * @param ck: element's third index in the third cube.
 */
template <int N, int D, size_t NM, size_t NMH>
inline void
Cube_Product_Checker<N, D, NM, NMH>::
decode_indices_from_index(int index,
                          int& ai, int& aj, int& ak,
                          int& bi, int& bj, int& bk,
                          int& ci, int& cj, int& ck) const
{
    ck = index % length;
    index /= length;
    cj = index % length;
    index /= length;
    ci = index % length;
    index /= length;

    bk = index % length;
    index /= length;
    bj = index % length;
    index /= length;
    bi = index % length;
    index /= length;

    ak = index % length;
    index /= length;
    aj = index % length;
    index /= length;
    ai = index;
    return;
}

//=============================================================================

/**
 * Get element's index in a matrix (2-dimensional case).
 *
 * @param i: element's first index in the cube;
 * @param j: element's second index in the cube;
 *
 * @return element's index in the cube.
 */
template <int N, int D, size_t NM, size_t NMH>
inline int
Cube_Product_Checker<N, D, NM, NMH>::
get_element_index(int i, int j) const
{
    return ((i*length) + j);
}

/**
 * Get element's index in a cube (3-dimensional case).
 *
 * @param i: element's first index in the cube;
 * @param j: element's second index in the cube;
 * @param k: element's third index in the cube;
 *
 * @return element's index in the cube.
 */
template <int N, int D, size_t NM, size_t NMH>
inline int
Cube_Product_Checker<N, D, NM, NMH>::
get_element_index(int i, int j, int k) const
{
    return (((i*length) + j)*length + k);
}

//=============================================================================

/**
 * Get index of a multiplication vector by its sum indices.
 *
 * This index is a unique number of the vector
 * in the set of all non-zero multiplication vectors.
 *
 * @param i: vector's first index;
 * @param j: vector's second index;
 *
 * @return multiplication vector index.
 */
template <int N, int D, size_t NM, size_t NMH>
inline int
Cube_Product_Checker<N, D, NM, NMH>::
get_m_index(int i, int j) const
{
    return ((i*m_length) + j);
}

/**
 * Get index of a multiplication vector by its sum indices.
 *
 * This index is a unique number of the vector
 * in the set of all non-zero multiplication vectors.
 *
 * @param i: vector's first index;
 * @param j: vector's second index;
 * @param k: vector's third index;
 *
 * @return multiplication vector index.
 */
template <int N, int D, size_t NM, size_t NMH>
inline int
Cube_Product_Checker<N, D, NM, NMH>::
get_m_index(int i, int j, int k) const
{
    return (((i*m_length) + j)*m_length + k);
}

//=============================================================================

/**
 * Get coefficients of cube elements by number of multiplication vector.
 *
 * @param multiplication vector index.
 *
 * @return
 */
template <int N, int D, size_t NM, size_t NMH>
inline vector<int>
Cube_Product_Checker<N, D, NM, NMH>::
decode_m_index(int index) const
{
    int last, previous;
    last = index % m_length;
    index /= m_length;
    previous = index % m_length;
    vector<int> result;
    if (dimension == 3) {
        int first = index / m_length;
        result.push_back(first);
    }
    result.push_back(previous);
    result.push_back(last);
    return result;
}

//=============================================================================

/**
 * Calculate product result vectors (2x2 case).
 */
template <>
void
Cube_Product_Checker<2, 2, 16, 4>::
calculate_r_vectors()
{
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j) {
            int index = get_element_index(i, j);
            r_vectors[index].v.reset();
            for (int l = 0; l < length; ++l) {
                r_vectors[index].v[get_vector_index(i, l, l, j)] = 1;
            }
        }
    }
}

/**
 * Calculate product result vectors (3x3 case).
 */
template <>
void
Cube_Product_Checker<3, 2, 81, 9>::
calculate_r_vectors()
{
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j) {
            int index = get_element_index(i, j);
            r_vectors[index].v.reset();
            for (int l = 0; l < length; ++l) {
                r_vectors[index].v[get_vector_index(i, l, l, j)] = 1;
            }
        }
    }
}

/**
 * Calculate product result vectors (2x2x2 case).
 */
template <>
void
Cube_Product_Checker<2, 3, 512, 8>::
calculate_r_vectors()
{
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j) {
            for (int k = 0; k < length; ++k) {
                int index = get_element_index(i, j, k);
                r_vectors[index].v.reset();
                for (int l = 0; l < length; ++l) {
                    r_vectors[index].v[get_vector_index(i, j, l, i, l, k, l, j, k)] = 1;
                }
            }
        }
    }
}

//=============================================================================

/**
 * Calculate non-zero multiplication vectors (2x2 case).
 */
template <>
void
Cube_Product_Checker<2, 2, 16, 4>::
calculate_m_vectors()
{
    for (int i = 1; i < power(2,element_count); ++i) {
        for (int j = 1; j < power(2,element_count); ++j) {
            Multiplication_Part_Vector av(i);
            Multiplication_Part_Vector bv(j);
            int index = get_m_index(i-1, j-1);
            m_vectors[index].v.reset();
            for (int k = 0; k < element_count; ++k) {
                for (int l = 0; l < element_count; ++l) {
                    if (av[k] && bv[l]) {
                        m_vectors[index].v.set(get_vector_index(k, l));
                    }
                }
            }
        }
    }
}

/**
 * Calculate non-zero multiplication vectors (3x3 case).
 */
template <>
void
Cube_Product_Checker<3, 2, 81, 9>::
calculate_m_vectors()
{
    for (int i = 1; i < power(2,element_count); ++i) {
        for (int j = 1; j < power(2,element_count); ++j) {
            Multiplication_Part_Vector av(i);
            Multiplication_Part_Vector bv(j);
            int index = get_m_index(i-1, j-1);
            m_vectors[index].v.reset();
            for (int k = 0; k < element_count; ++k) {
                for (int l = 0; l < element_count; ++l) {
                    if (av[k] && bv[l]) {
                        m_vectors[index].v.set(get_vector_index(k, l));
                    }
                }
            }
        }
    }
}

/**
 * Calculate non-zero multiplication vectors (3-d case).
 */
template <>
void
Cube_Product_Checker<2, 3, 512, 8>::
calculate_m_vectors()
{
    for (int i = 1; i < power(2,element_count); ++i) {
        for (int j = 1; j < power(2,element_count); ++j) {
            for (int k = 1; k < power(2,element_count); ++k) {
                Multiplication_Part_Vector av(i);
                Multiplication_Part_Vector bv(j);
                Multiplication_Part_Vector cv(k);
                int index = get_m_index(i-1, j-1, k-1);
                m_vectors[index].v.reset();
                for (int l = 0; l < element_count; ++l) {
                    for (int o = 0; o < element_count; ++o) {
                        for (int p = 0; p < element_count; ++p) {
                            if (av[l] & bv[o] & cv[p]) {
                                m_vectors[index].v.set(get_vector_index(l, o, p));
                            }
                        }
                    }
                }
            }
        }
    }
}

//=============================================================================

/**
 * Clear the current sets.
 */
template <int N, int D, size_t NM, size_t NMH>
void
Cube_Product_Checker<N, D, NM, NMH>::
clear_sets()
{
    n_vectors_indexes.clear();
    good_vectors_indexes.clear();
}

//=============================================================================

/**
 * Add a multiplication vector to the current set.
 *
 * @param index: index of the multiplication vector.
 */
template <int N, int D, size_t NM, size_t NMH>
void
Cube_Product_Checker<N, D, NM, NMH>::
add_vector_to_set(int index)
{
    n_vectors_indexes.insert(index);
}

//=============================================================================

/**
 * Clear statistics.
 */
template <int N, int D, size_t NM, size_t NMH>
void
Cube_Product_Checker<N, D, NM, NMH>::
clear_statistics()
{
    checked_sets_count = 0;
    iteration_count = 0;
    best_result = 0;
    lin_dependent_sets = 0;
#ifdef OUTPUT_STATISTICS
    raw_solution_count = 0;
    solutions.clear();
    solution_distribution.clear();
#endif // OUTPUT_STATISTICS
}

//=============================================================================

/**
 * Make a random candidate.
 */
template <int N, int D, size_t NM, size_t NMH>
void
Cube_Product_Checker<N, D, NM, NMH>::
make_random_candidate()
{
    clear_sets();
    for (int i = 0; i < (f_count-element_count); ++i) {
        int cc = rnd.next();
        while (n_vectors_indexes.count(cc) > 0) {
            cc = rnd.next();
        }
        add_vector_to_set(cc);
    }
}

//=============================================================================

/**
 * Make a random candidate.
 */
template <int N, int D, size_t NM, size_t NMH>
void
Cube_Product_Checker<N, D, NM, NMH>::
make_candidate_from_space()
{
    clear_sets();
    n_vectors_indexes = (*candidate_space.begin());
}

//=============================================================================

/**
 * If current candidate is in cache, than get info from cache.
 *
 * @return true if the current candidate is in cache.
 */
template <int N, int D, size_t NM, size_t NMH>
bool
Cube_Product_Checker<N, D, NM, NMH>::
check_cache()
{
    map<Candidate, int>::const_iterator it;
    if ((it = candidate_cache.find(n_vectors_indexes)) != candidate_cache.end()) {
        best_result = it->second;
        ++cache_hits;
        return true;
    } else {
        return false;
    }
}

//=============================================================================

/**
 * Update the cache with new result.
 */
template <int N, int D, size_t NM, size_t NMH>
void
Cube_Product_Checker<N, D, NM, NMH>::
update_cache()
{
    if (candidate_cache.size() >= cache_limit) {
        candidate_cache.clear();
        //cout << "------------->> Cache cleared" << endl;
    }
    candidate_cache[n_vectors_indexes] = best_result;
}

//=============================================================================

/**
 * Check the current set for being a solution.
 *
 * @return true if the current set is a solution.
 */
template <int N, int D, size_t NM, size_t NMH>
bool
Cube_Product_Checker<N, D, NM, NMH>::
check_vectors_for_goodness()
{
#ifdef VERY_DETAILED_OUTPUT
    cout << "<" << thread_number << "> ";
    cout << "Checking of vectors { ";
    for (int i: n_vectors_indexes) {
        cout << i << " ";
    }
    cout << "} has started..." << endl;
    tw.watch();
#endif // VERY_DETAILED_OUTPUT
    ++checked_sets_count;
#ifdef USE_CACHE
    if (check_cache()) {
        if (best_result == f_count) {
            return true;
        } else {
            return false;
        }
    }
#endif // USE_CACHE
    vector<mm_vector_with_properties<NM>> nvwp; /// set of span vectors for SLAE
    vector<mm_vector_with_properties<NM>> gvwp; /// set of good vectors for SLAE
    Vectors_Presolve_Data<NM> v; /// vector presolve data for span vectors
    Gauss_WP_Presolve_Data<NM> pwp(0); /// presolve data for SLAE for span vectors
    Gauss_WP_Presolve_Data<NM> pwpg(0); /// presolve data for SLAE for good vectors

    for (set<int>::iterator i = n_vectors_indexes.begin(); i != n_vectors_indexes.end(); ++i) {
        nvwp.push_back(m_vectors[*i]);
        v.add_vector(m_vectors[*i].v);
    }
    for (int i = 0; i < element_count; ++i) {
        nvwp.push_back(r_vectors[i]);
        v.add_vector(r_vectors[i].v);
    }
    if (!gauss_wp_presolve(nvwp, pwp)) { // vectors are linearly dependent
        ++lin_dependent_sets;
#ifdef USE_CACHE
        update_cache();
#endif // USE_CACHE
        return false;
    }
    gauss_wp_presolve(gvwp, pwpg);
    v.presolve();
    for (int i = 0; i < m_count; ++i) {
        //if (!v.check(m_vectors[i].v)) { // discard vector by checking bits
        //    ++bit_check_hits;
        //    continue;
        //}
#ifdef OUTPUT_STATISTICS
        ++gaussian_eliminations;
#endif // OUTPUT_STATISTICS
        if (gauss_wp_solve(pwp, m_vectors[i])) { // if vector is in the current span
#ifdef OUTPUT_STATISTICS
            ++gaussian_eliminations;
#endif // OUTPUT_STATISTICS
            if (!gauss_wp_solve(pwpg, m_vectors[i])) { // if not linearly dependent
                good_vectors_indexes.insert(i);
                gvwp.push_back(m_vectors[i]);
                gauss_wp_presolve(gvwp, pwpg);
            }
            if (good_vectors_indexes.size() >= f_count) // we have found solution
                break;
        }
    }
#ifdef VERY_DETAILED_OUTPUT
    cout << "  [" << tw.watch() << " s] Done." << endl;
#endif // VERY_DETAILED_OUTPUT
    best_result = good_vectors_indexes.size();
    if (best_result > (f_count-element_count)) {
        if (best_result >= top_best_result) {
            if (best_result > top_best_result) {
                top_best_result = best_result;
                top_best_solutions.clear();
#ifdef SAVE_BEST_RESULTS_TO_FILE
                ofstream fout(to_string(thread_number)+string("bestsofar.txt"), ios_base::app);
                fout << top_best_result << " vectors" << endl;
                fout.close();
#endif // SAVE_BEST_RESULTS_TO_FILE
            }
            top_best_solutions.insert(n_vectors_indexes);
#ifdef SAVE_BEST_RESULTS_TO_FILE
            ofstream fout(to_string(thread_number)+string("bestsofar.txt"), ios_base::app);
            for (set<int>::iterator cc = n_vectors_indexes.begin(); cc != n_vectors_indexes.end(); ++cc) {
                fout << *cc << " ";
            }
            fout << endl;
            fout.close();
#endif // SAVE_BEST_RESULTS_TO_FILE
        }
    }
#ifdef USE_CACHE
        update_cache();
#endif // USE_CACHE
    if (good_vectors_indexes.size() >= f_count) { // there is a solution
#ifdef VERBOSE_OUTPUT
        cout << "  Good vectors have been found: { ";
        for (set<int>::iterator cc = good_vectors_indexes.begin(); cc != good_vectors_indexes.end(); ++cc) {
            cout << *cc << " ";
        }
        cout << "}" << endl;
#endif // VERBOSE_OUTPUT
#ifdef OUTPUT_SOLUTIONS_TO_FILE
        ofstream mvfile("good.txt");
        mvfile << "Found good vectors!!!\n";
        for (set<int>::iterator cc = good_vectors_indexes.begin(); cc != good_vectors_indexes.end(); ++cc) {
            mvfile << *cc << " ";
        }
        mvfile << "\n";
#endif // OUTPUT_SOLUTIONS_TO_FILE
        return true;
    }
    return false;
}

//=============================================================================

/**
 * Check all solution space for solutions (2x2 case).
 *
 * @return true if at least one solution was found.
 */
template <>
bool
Cube_Product_Checker<2, 2, 16, 4>::
check_for_good_vectors()
{
    clear_statistics();
    for (int c1 = 0; c1 < m_count-2; ++c1) {
        //cout << c1 << endl;
        for (int c2 = c1+1; c2 < m_count-1; ++c2) {
            for (int c3 = c2+1; c3 < m_count; ++c3) {
                if (*stop_signal)
                    return false;
                clear_sets();
                add_vector_to_set(c1);
                add_vector_to_set(c2);
                add_vector_to_set(c3);
                if (check_vectors_for_goodness()) {
#ifdef OUTPUT_STATISTICS
                    solutions.insert(good_vectors_indexes);
                    ++solution_distribution[good_vectors_indexes];
                    ++raw_solution_count;
#endif // OUTPUT_STATISTICS
                }
            }
        }
    }
    return true;
}

/**
 * Check all solution space for solutions (2x2x2 case).
 *
 * @return true if a solution was found.
 */
template <>
bool
Cube_Product_Checker<2, 3, 512, 8>::
check_for_good_vectors()
{
    clear_statistics();
    for (int c1 = 0; c1 < m_count-6; ++c1)
        for (int c2 = c1+1; c2 < m_count-5; ++c2)
            for (int c3 = c2+1; c3 < m_count-4; ++c3)
                for (int c4 = c3+1; c4 < m_count-3; ++c4)
                    for (int c5 = c4+1; c5 < m_count-2; ++c5)
                        for (int c6 = c5+1; c6 < m_count-1; ++c6)
                            for (int c7 = c6+1; c7 < m_count; ++c7) {
                                if (*stop_signal)
                                    return false;
                                clear_sets();
                                add_vector_to_set(c1);
                                add_vector_to_set(c2);
                                add_vector_to_set(c3);
                                add_vector_to_set(c4);
                                add_vector_to_set(c5);
                                add_vector_to_set(c6);
                                add_vector_to_set(c7);
                                if (check_vectors_for_goodness())
                                    return true;
                            }
    return false;
}

//=============================================================================

/**
 * Check solutions by randomly picking vector sets.
 *
 * @return true if a solution was found.
 */
template <int N, int D, size_t NM, size_t NMH>
bool
Cube_Product_Checker<N, D, NM, NMH>::
check_for_good_vectors_randomized()
{
    clear_statistics();
    while (true) {
        if (*stop_signal)
            return false;
        make_random_candidate();
        if (check_vectors_for_goodness())
            return true;
    }
    return false;
}

//=============================================================================

/**
 * Make neighbourhood of the current candidate.
 */
template <int N, int D, size_t NM, size_t NMH>
void
Cube_Product_Checker<N, D, NM, NMH>::
start_neighbourhood()
{
    //cout << "Start making neighbours" << endl;
    neighbours.clear();
    // all possible one-coefficient changes in one of the vectors
    for (int i: n_vectors_indexes) {
        //cout << "  start with vector " << i << endl;
        Candidate c_set = Candidate(n_vectors_indexes);
        c_set.erase(i);
        vector<int> coef = decode_m_index(i);
        vector<Multiplication_Part_Vector> cur_coef;
        for (int j = 0; j < dimension; ++j) {
            cur_coef.push_back(Multiplication_Part_Vector(coef[j]+1));
        }
#ifdef OUTPUT_STATISTICS
        int value = cur_coef[0].to_ulong()-1;
        for (int o = 1; o < dimension; ++o) {
            value = value*m_length + cur_coef[o].to_ulong()-1;
        }
        if (value != i) {
            cout << "Hardcore error here! " << value << endl;
        }
#endif // OUTPUT_STATISTICS
        for (int j = 0; j < dimension; ++j) {
            //cout << "    start with cube " << j << " with total bits " << cur_coef[j] << endl;
            for (int k = 0; k < cur_coef[j].size(); ++k) {
                //cout << "      start with bit " << k << endl;
                cur_coef[j].flip(k);
                if (cur_coef[j].count() == 0) {
                    cur_coef[j].flip(k);
                    continue;
                }
                Candidate new_set = Candidate(c_set);
                int value = cur_coef[0].to_ulong()-1;
                for (int o = 1; o < dimension; ++o) {
                    value = value*m_length + cur_coef[o].to_ulong()-1;
                }
#ifdef OUTPUT_STATISTICS
                if ((value < 0) || (value >= m_count)) {
                    cout << "--->  got value " << value << endl
                              << "--->  after " << i << endl
                              << "--->  on j = " << j << " and k = " << k << endl
                              << "--->  with cur_coef[j] = " << cur_coef[j] << " and it's long as " << cur_coef[j].to_ulong() << endl;
                }
#endif // OUTPUT_STATISTICS
                new_set.insert(value);
                if (new_set.size() == n_vectors_indexes.size()) {
                    neighbours.insert(new_set);
                }
                cur_coef[j].flip(k);
            }
        }
    }
    //cout << "Made " << neighbours.size() << " neigbours." << endl;
}

//=============================================================================

/**
 * Check solutions by doing local search.
 *
 * @return true if a solution was found.
 */
template <int N, int D, size_t NM, size_t NMH>
bool
Cube_Product_Checker<N, D, NM, NMH>::
solve_hill_climbing(int local_max_limit, bool use_space)
{
    Random r;
    clear_statistics();
    clear_sets();
    if (use_space) {
        make_candidate_from_space();
    } else {
        make_random_candidate();
    }
    check_vectors_for_goodness();
    local_max_iterations = 0;
    local_iterations = 0;
    restarts = 0;
    while (true) {
        if (local_max_iterations > local_max_limit) {
#ifdef VERBOSE_OUTPUT
            cout << omp_get_thread_num() << "] Doing restart! local iterations = " << local_iterations << " / " << iteration_count << endl;
#endif // VERBOSE_OUTPUT
            clear_sets();
            if (use_space) {
                make_candidate_from_space();
            } else {
                make_random_candidate();
            }
            check_vectors_for_goodness();
            local_max_iterations = 0;
            local_iterations = 0;
            ++restarts;
        }
#ifdef VERBOSE_OUTPUT
        cout << omp_get_thread_num() << "] Iteration " << local_iterations << ": best is " << best_result << endl;
#endif // VERBOSE_OUTPUT
        if (best_result >= f_count)
            break;
        int next_best_result = 0;
        int old_best_result = best_result;
        vector<Candidate> best_candidates;
        start_neighbourhood();
        for (const Candidate& c: neighbours) {
            if (*stop_signal)
                return false;
            n_vectors_indexes = c;
            good_vectors_indexes.clear();
            check_vectors_for_goodness();
            if (best_result > next_best_result) {
                next_best_result = best_result;
                best_candidates.clear();
                best_candidates.push_back(n_vectors_indexes);
            } else if (best_result == next_best_result) {
                best_candidates.push_back(n_vectors_indexes);
            }
        }
        r.init(0, best_candidates.size()-1);
        n_vectors_indexes = best_candidates[r.next()];
        best_result = next_best_result;
        if (next_best_result > old_best_result) {
            local_max_iterations = 0;
        } else {
            ++local_max_iterations;
        }
        ++iteration_count;
        ++local_iterations;
#ifdef VERBOSE_OUTPUT
        output_current_state();
#endif // VERBOSE_OUTPUT
        //cout << "Iteration " << iteration_count << ": best is " << best_result << endl;
    }
    good_vectors_indexes.clear();
    check_vectors_for_goodness();
    if (best_result >= f_count)
        return true;
    else {
        cout << "Error!" << endl;
        return false;
    }
}

//=============================================================================

/**
 * Output current state variables of the checker.
 */
template <int N, int D, size_t NM, size_t NMH>
void
Cube_Product_Checker<N, D, NM, NMH>::
output_current_state() const
{
    cout << endl;
    cout << "===================================================" << endl;
    cout << "Current state of thread <" << thread_number << "> is:" << endl;
    cout << checked_sets_count << " checked sets" << endl;
    cout << restarts << " restarts" << endl;
    cout << iteration_count << " iterations" << endl;
    cout << lin_dependent_sets << " linearly dependent sets hits" << endl;
    cout << bit_check_hits << " bit check hits" << endl;
    cout << gaussian_eliminations << " gaussian eliminations" << endl;
    cout << cache_hits << " cache hits" << endl;
    cout << "===================================================" << endl;
}

//=============================================================================

/**
 * Write m_vectors into file.
 *
 * @param filename: filename.
 */
template <int N, int D, size_t NM, size_t NMH>
void
Cube_Product_Checker<N, D, NM, NMH>::
write_m_vectors(const string& filename)
{
    ofstream fout(filename);
    for (int i = 0; i < m_count; ++i) {
        fout << m_vectors[i].v.to_string() << m_vectors[i].r.to_string();
    }
    fout.close();
}

//=============================================================================

/**
 * Read m_vectors from file.
 *
 * @param filename: filename.
 */
template <int N, int D, size_t NM, size_t NMH>
void
Cube_Product_Checker<N, D, NM, NMH>::
read_m_vectors(const string& filename)
{
    ifstream fin(filename);
    for (int i = 0; i < m_count; ++i) {
        fin >> m_vectors[i].v >> m_vectors[i].r;
    }
    fin.close();
}

//=============================================================================

/**
 * Read candidate space from file.
 *
 * @param filename: filename.
 */
template <int N, int D, size_t NM, size_t NMH>
void
Cube_Product_Checker<N, D, NM, NMH>::
read_candidate_space(const string& filename)
{
    ifstream fin(filename);
    set<int> candidate;
    for (int i = 0; i < f_count-element_count; ++i) {
        int v;
        fin >> v;
        candidate.insert(v);
    }
    candidate_space.insert(candidate);
    fin.close();
}

//=============================================================================

/**
 * Print binary value of a vector.
 *
 * @param v: the vector.
 */
template <int N, int D, size_t NM, size_t NMH>
void
Cube_Product_Checker<N, D, NM, NMH>::
output_vector(Multiplication_Vector v) const
{
    cout << v;
}

//=============================================================================

/**
 * Print value of a vector in letters (2-d case).
 *
 * @param v: the vector.
 */
template <>
void
Cube_Product_Checker<2, 2, 16, 4>::
output_vector_text(Multiplication_Vector v) const
{
    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i]) {
            int ai, aj, bi, bj;
            decode_indices_from_index(i, ai, aj, bi, bj);
            cout << "A" << ai+1 << aj+1 << "B" << bi+1 << bj+1 << " ";
        }
    }
}

//=============================================================================

/**
 * Generate and save to a file a sequence of random sets.
 *
 * @param size: number of sets;
 * @param filename: filename.
 */
template <int N, int D, size_t NM, size_t NMH>
void
Cube_Product_Checker<N, D, NM, NMH>::
save_random_samples(int size, const char* filename) const
{
    ofstream fout(filename);
    Random rnd(0, m_count-1);
    fout << size << "\n";
    for (int i = 0; i < size; ++i) {
        clear_sets();
        for (int i = 0; i < (f_count-element_count); ++i) {
            int cc = rnd.next();
            while (n_vectors_indexes.count(cc) > 0) {
                cc = rnd.next();
            }
            add_vector_to_set(cc);
            fout << cc << " ";
        }
        fout << "\n";
    }
    fout.close();
}

//=============================================================================

/**
 * Read a set sequence from the file and check it.
 *
 * @param filenamein: filename to read from;
 * @param filenameout: filename for output.
 */
template <int N, int D, size_t NM, size_t NMH>
void Cube_Product_Checker<N, D, NM, NMH>::
read_samples_and_check(const char* filenamein, const char* filenameout) const
{
    ifstream fin(filenamein);
    ofstream fout(filenameout, ios_base::app);
    fout << "\n ################################### \n";
    int size;
    fin >> size;
    Timewatch timer;
    double time = 0.0;
    int gv = 0;
    for (int i = 0; i < size; ++i) {
        if (N*D > 5) {
            cout << "working on case " << i << "\n";
        }
        clear_sets();
        for (int i = 0; i < (f_count-element_count); ++i) {
            int cc;
            fin >> cc;
            add_vector_to_set(cc);
        }
        timer.watch();
        if (check_vectors_for_goodness()) {
            ++gv;
        }
        double curtime = timer.watch();
        time += curtime;
        if (N*D > 5) {
            fout << "\t" << i << ": " << curtime << " s\n";
        }
    }
    cout << "Found " << gv << " solutions\n";
    fout << "Total time: " << time << " s (" << (time/size) << " s avg) (found " << gv << " good vectors)\n";
    fout.close();
}

//=============================================================================

/**
 * Calculate statistics and save to a file.
 *
 * @param filename: filename to save statistics.
 */
template <int N, int D, size_t NM, size_t NMH>
void
Cube_Product_Checker<N, D, NM, NMH>::
save_results(const char* filename) {
    ofstream fout(filename);
    fout << "Results: " << solutions.size() << " different solutions were found." << endl;
    fout << raw_solution_count << " vector sets were successful." << endl;
    map<int, int> vector_distribution = map<int, int>();
    vector<map<int, int>> bit_distributions;
    for (set<set<int>>::iterator i = solutions.begin(); i != solutions.end(); ++i) {
        fout << "[ ";
        map<int, int> bd;
        for (set<int>::iterator j = (*i).begin(); j != (*i).end(); ++j) {
            fout << (*j) << " ";
            ++vector_distribution[*j];
            for (size_t k = 0; k < m_vectors[*j].v.size(); ++k) {
                if (m_vectors[*j].v.test(k)) {
                    ++bd[k];
                }
            }
        }
        bit_distributions.push_back(bd);
        fout << "] : " << solution_distribution[(*i)] << endl;
    }
    fout << endl << "Used vectors are:" << endl;
    for (const auto& e: vector_distribution) {
        fout << e.first << "\t" << m_vectors[e.first].v << " : " << e.second << endl;
    }
    fout << endl << "Bit distributions per solution are:" << endl;
    for (auto& e: bit_distributions) {
        for (int i = 0; i < NM; ++i) {
            fout << e[i];
        }
        fout << endl;
    }
    fout.close();
}

//=============================================================================

/**
 * Check if a given set of vectors is a valid solution.
 *
 * @param s:  set of result vectors;
 *
 * @param sp: solution properties;
 *
 * @return if it is a solution.
 */
template <int N, int D, size_t NM, size_t NMH>
bool
Cube_Product_Checker<N, D, NM, NMH>::
check_solution(set<int> s, Solution_Properties& sp) {
    sp.coefficients.clear();
    sp.multiplication_vectors.clear();
    sp.multiplication_vectors = s;
    sp.operation_count = 0;
    vector<mm_bitset<NM>> vectors;
    for (auto i: s) {
        vectors.push_back(m_vectors[i].v);
        vector<int> ops = sum_operations_cube(i);
        sp.operation_count += accumulate(ops.begin(), ops.end(), -ops.size());
    }
    for (int i = 0; i < element_count; ++i) {
        if (!gauss_solve(vectors, r_vectors[i].v)) {
            return false;
        }
        boost::dynamic_bitset<> x = binary_solve_result(vectors, r_vectors[i].v);
        sp.coefficients.push_back(x);
        sp.operation_count += x.count() - 1;
    }
    return true;
}

//=============================================================================

/**
 * Return number of summation operations inside cubes for a given index of multiplication vector.
 *
 * @param index: multiplication vector index;
 *
 * @return list of operation count for each cube.
 */
template <int N, int D, size_t NM, size_t NMH>
vector<int>
Cube_Product_Checker<N, D, NM, NMH>::
sum_operations_cube(int index)
{
    vector<int> result;
    if (dimension == 2) {
        int i, j;
        j = index % m_length + 1;
        i = index / m_length + 1;
        result.push_back(popcount(i));
        result.push_back(popcount(j));
    } else if (dimension == 3) {
        int i, j, k;
        k = index % m_length + 1;
        index /= m_length;
        j = index % m_length + 1;
        i = index / m_length + 1;
        result.push_back(popcount(i));
        result.push_back(popcount(j));
        result.push_back(popcount(k));
    }
    return result;
}

//=============================================================================

/**
 * Output solution properties to a file in readable text form.
 *
 * @param sp: solution properties object;
 * @param filename: filename.
 */
template <int N, int D, size_t NM, size_t NMH>
void
Cube_Product_Checker<N, D, NM, NMH>::
save_solution_properties(const Solution_Properties& sp, const char* filename)
{
    ofstream out(filename);
    //----- header
    out << "Strassen solution for ";
    for (int i = 0; i < dimension-1; ++i) {
        out << length << "×";
    }
    out << length << " cube product with " << f_count << " multiplications." << endl << endl;
    //----- multiplication vectors numbers
    out << "Multiplication numbers: [ ";
    for (auto i: sp.multiplication_vectors) {
        out << i << " ";
    }
    out << "]" << endl;
    //----- multiplications
    int c = 1;
    int m_summations = 0;
    for (auto i = sp.multiplication_vectors.begin(); i != sp.multiplication_vectors.end(); ++i, ++c) {
        out << "M" << c << " = ";
        vector<int> ops = sum_operations_cube(*i);
        m_summations += accumulate(ops.begin(), ops.end(), -ops.size());
        vector<string> sm, smc;
        set<string> sa, sb, sc;
        Multiplication_Vector& mv = m_vectors[*i].v;
        if (dimension == 2) {
            for (size_t j = 0; j < mv.size(); ++j) {
                if (mv[j]) {
                    int ai, aj, bi, bj;
                    decode_indices_from_index(j, ai, aj, bi, bj);
                    sm.push_back("A" + to_string(ai+1) + to_string(aj+1) +
                                 "B" + to_string(bi+1) + to_string(bj+1));
                    sa.insert("A" + to_string(ai+1) + to_string(aj+1));
                    sb.insert("B" + to_string(bi+1) + to_string(bj+1));
                }
            }
            out << "(" << boost::algorithm::join(sa, " + ") + ")×";
            out << "(" << boost::algorithm::join(sb, " + ") + ")";
        } else if (dimension == 3) {
            // TODO
        }
        out << " = ";
        out << boost::algorithm::join(sm, " + ") << endl;
    }
    //----- result elements
    int r_summations = 0;
    for (int i = 0; i < element_count; ++i) {
        vector<string> sr;
        boost::dynamic_bitset<> x = sp.coefficients[i];
        r_summations += x.count() - 1;
        if (dimension == 2) {
            out << "C" << i/length+1 << i%length+1 << " = ";
        } else if (dimension == 3) {
            // TODO
        }
        for (boost::dynamic_bitset<>::size_type j = 0; j < x.size(); ++j) {
            if (x.test(j)) {
                sr.push_back("M" + to_string(j+1));
            }
        }
        out << boost::algorithm::join(sr, " + ") << endl;
    }
    //----- summation count
    out << "Summation count for multiplications = " << m_summations << endl;
    out << "Summation count for result elements = " << r_summations << endl;
    out << "Total number of summations used = " << sp.operation_count << endl;
    //----- end
    out.close();
}

//=============================================================================

#endif // MMCHECKER_HPP_INCLUDED
