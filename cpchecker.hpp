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

#include "slae.hpp"
#include "utils.hpp"

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
    int length; /// size of cube (N)
    int dimension; /// dimension of cube (D)
    int element_count; /// number of elements in cube (N^D)
    mm_vector_with_properties_options vector_options; /// options for vector properties
    mm_vector_with_properties<NM>* r_vectors; /// result vectors
    mm_vector_with_properties<NM>* m_vectors; /// multiplication vectors
    int m_count; /// number of non-zero multiplication vectors = (2^(N^D)-1)^D
    int m_length; /// number non-zero element sums = 2^(N^D)-1
    int f_count; /// number of vectors to choose for a basis = N^(D+1)-1
    std::set<int> good_vectors_indexes; /// set of indexes for vectors that are in the current span
    std::set<int> n_vectors_indexes; /// set of indexes of current chosen vectors
    std::set<std::set<int>> solutions; /// set of found solutions
    Timewatch tw; /// timer that is used for getting calculation time

    //=============--- constructors and destructors
    Cube_Product_Checker();
    ~Cube_Product_Checker();

    //=============--- index operations
    /// return linear index in a bit vector by its indexes in matrices
    int get_vector_index(int ai, int aj, int bi, int bj); /// for 2-dimensional case
    int get_vector_index(int ai, int aj, int ak, int bi, int bj, int bk, int ci, int cj, int ck); /// for 3-d case
    /// return linear index in a bit vector by combined indexes in matrices
    int get_vector_index(int a, int b); /// for 2-d case
    int get_vector_index(int a, int b, int c); /// for 3-d case
    /// return indices from index
    void decode_indices_from_index(int index, int& ai, int& aj, int& bi, int& bj); /// for 2-d case
    void decode_indices_from_index(int index, int& ai, int& aj, int& ak, int& bi, int& bj, int& bk, int& ci, int& cj, int& ck); /// for 3-d case
    /// return linear index of an element in a matrix
    int get_element_index(int i, int j); /// for 2-d case
    int get_element_index(int i, int j, int k); /// for 3-d case
    /// return index of multiplication vector in the set
    int get_m_index(int i, int j); /// for 2-d case
    int get_m_index(int i, int j, int k); /// for 3-d case

    //=============--- Initial calculations
    void init(); /// calculate all properties
    void calculate_r_vectors(); /// write result vectors to array
    void calculate_m_vectors(); /// write multiplication vectors to array

    //=============--- checking routines
    void clear_sets(); /// clear current sets of vectors
    void add_vector_to_set(int index); /// add vector to current sets
    bool check_vectors_for_goodness(); /// check current set of vectors

    //=============--- searching for solution
    bool check_for_good_vectors(); /// check all solution space
    bool check_for_good_vectors_randomized(); /// do random search

    //=============--- utilities
    void output_vector(Multiplication_Vector v); /// output vector to screen
    void output_vector_text(Multiplication_Vector v); /// output vector to screen in letters
    void save_random_samples(int size, const char* filename); /// save random sets to a file (for testing later)
    void read_samples_and_check(const char* filename, const char* filenameout); /// check sets from a file
    bool check_vectors_for_goodness_a1(); /// just normal gauss everywhere
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
Cube_Product_Checker()
{
    length = N;
    dimension = D;
    element_count = power(length, dimension);
    f_count = power(length,dimension+1)-1;
    r_vectors = new mm_vector_with_properties<NM>[element_count];
    m_length = power(2,element_count)-1;
    m_count = power(m_length,dimension);
    m_vectors = new mm_vector_with_properties<NM>[m_count];
    mm_vector_with_properties<NM>::make_options(vector_options);
#ifdef VERBOSE_OUTPUT
    std::cout << "Cube Product Checker has been created." << std::endl;
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
    delete [] r_vectors;
    delete [] m_vectors;
}

//=============================================================================

/**
 * Initialize the object with all necessary properties.
 */
template <int N, int D, size_t NM, size_t NMH>
void
Cube_Product_Checker<N, D, NM, NMH>::
init()
{
#ifdef VERBOSE_OUTPUT
    tw.watch();
#endif // VERBOSE_OUTPUT
    calculate_r_vectors();
    for (int i = 0; i < element_count; ++i) {
        r_vectors[i].calculate_properties(vector_options);
    }
#ifdef VERBOSE_OUTPUT
    std::cout << "[" << tw.watch() << " s] Result vectors calculated." << std::endl;
#endif // VERBOSE_OUTPUT
    calculate_m_vectors();
    for (int i = 0; i < m_count; ++i) {
        m_vectors[i].calculate_properties(vector_options);
    }
#ifdef VERBOSE_OUTPUT
    std::cout << "[" << tw.watch() << " s] Multiplication vectors calculated." << std::endl;
#endif // VERBOSE_OUTPUT
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
                 int bi, int bj)
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
                 int ci, int cj, int ck)
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
get_vector_index(int a, int b)
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
get_vector_index(int a, int b, int c)
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
                          int& bi, int& bj)
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
                          int& ci, int& cj, int& ck)
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
get_element_index(int i, int j)
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
get_element_index(int i, int j, int k)
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
get_m_index(int i, int j)
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
get_m_index(int i, int j, int k)
{
    return (((i*m_length) + j)*m_length + k);
}

//=============================================================================

/**
 * Calculate product result vectors (2-d case).
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
 * Calculate product result vectors (3-d case).
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
 * Calculate non-zero multiplication vectors (2-d case).
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
    std::cout << "Checking of vectors { ";
    for (int i: n_vectors_indexes) {
        std::cout << i << " ";
    }
    std::cout << "} has started..." << std::endl;
    tw.watch();
#endif // VERY_DETAILED_OUTPUT
    std::vector<mm_vector_with_properties<NM>> nvwp; /// set of span vectors for SLAE
    std::vector<mm_vector_with_properties<NM>> gvwp; /// set of good vectors for SLAE
    Vectors_Presolve_Data<NM> v; /// vector presolve data for span vectors
    Gauss_WP_Presolve_Data<NM> pwp; /// presolve data for SLAE for span vectors
    Gauss_WP_Presolve_Data<NM> pwpg; /// presolve data for SLAE for good vectors

    for (std::set<int>::iterator i = n_vectors_indexes.begin(); i != n_vectors_indexes.end(); ++i) {
        nvwp.push_back(m_vectors[*i]);
        v.add_vector(m_vectors[*i].v);
    }
    for (int i = 0; i < element_count; ++i) {
        nvwp.push_back(r_vectors[i]);
        v.add_vector(r_vectors[i].v);
    }
    gauss_wp_presolve(nvwp, pwp);
    gauss_wp_presolve(gvwp, pwpg);

    for (int i = 0; i < m_count; ++i) {
        if (!v.check(m_vectors[i].v))
            continue; // discard vector by checking bits
        if (gauss_wp_solve(pwp, m_vectors[i])) { // if vector is in the current span
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
    std::cout << "  [" << tw.watch() << " s] Done." << std::endl;
#endif // VERY_DETAILED_OUTPUT
    if (good_vectors_indexes.size() >= f_count) { // there is a solution
#ifdef VERBOSE_OUTPUT
        std::cout << "  Good vectors have been found: { ";
        for (std::set<int>::iterator cc = good_vectors_indexes.begin(); cc != good_vectors_indexes.end(); ++cc) {
            std::cout << *cc << " ";
        }
        std::cout << "}" << std::endl;
#endif // VERBOSE_OUTPUT
#ifdef OUTPUT_SOLUTIONS_TO_FILE
        std::ofstream mvfile("good.txt");
        mvfile << "Found good vectors!!!\n";
        for (std::set<int>::iterator cc = good_vectors_indexes.begin(); cc != good_vectors_indexes.end(); ++cc) {
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
 * Check all solution space for solutions (2-d case).
 *
 * Saves solutions
 *
 * @return true if at least one solution was found.
 */
template <>
bool
Cube_Product_Checker<2, 2, 16, 4>::
check_for_good_vectors()
{
    solutions.clear();
    for (int c1 = 0; c1 < m_count-2; ++c1) {
        for (int c2 = c1+1; c2 < m_count-1; ++c2) {
            for (int c3 = c2+1; c3 < m_count; ++c3) {
                clear_sets();
                add_vector_to_set(c1);
                add_vector_to_set(c2);
                add_vector_to_set(c3);
                if (check_vectors_for_goodness()) {
                    solutions.insert(good_vectors_indexes);
                }
            }
        }
    }
#ifdef OUTPUT_STATISTICS
    std::ofstream fout("results2.txt");
    fout << "Results: " << results.size() << " found solutions.\n";
    for (std::set<std::set<int>>::iterator i = results.begin(); i != results.end(); ++i) {
        fout << "[ ";
        for (std::set<int>::iterator j = (*i).begin(); j != (*i).end(); ++j) {
            fout << (*j) << " ";
        }
        fout << "]\n";
    }
    fout.close();
#endif // OUTPUT_STATISTICS
    return (solutions.size() > 0);
}

/**
 * Check all solution space for solutions (3-d case).
 *
 * @return true if a solution was found.
 */
template <>
bool
Cube_Product_Checker<2, 3, 512, 8>::
check_for_good_vectors()
{
    solutions.clear();
    for (int c1 = 0; c1 < m_count-6; ++c1)
        for (int c2 = c1+1; c2 < m_count-5; ++c2)
            for (int c3 = c2+1; c3 < m_count-4; ++c3)
                for (int c4 = c3+1; c4 < m_count-3; ++c4)
                    for (int c5 = c4+1; c5 < m_count-2; ++c5)
                        for (int c6 = c5+1; c6 < m_count-1; ++c6)
                            for (int c7 = c6+1; c7 < m_count; ++c7) {
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
    Random rnd(0, m_count-1);
    while (true) {
        clear_sets();
        for (int i = 0; i < (f_count-element_count); ++i) {
            int cc = rnd.next();
            while (n_vectors_indexes.count(cc) > 0) {
                cc = rnd.next();
            }
            add_vector_to_set(cc);
        }
        if (check_vectors_for_goodness())
            return true;
    }
    return false;
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
output_vector(Multiplication_Vector v)
{
    std::cout << v;
}

//=============================================================================

/**
 * Print value of a vector in letters.
 *
 * @param v: the vector.
 */
template <int N, int D, size_t NM, size_t NMH>
void
Cube_Product_Checker<N, D, NM, NMH>::
output_vector_text(Multiplication_Vector v)
{
    for (int i = 0; i < v.size(); ++i) {
        if (v[i]) {
            int ai, aj, bi, bj;
            decode_indices_from_index(i, ai, aj, bi, bj);
            std::cout << "A" << ai+1 << aj+1 << "B" << bi+1 << bj+1 << " ";
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
save_random_samples(int size, const char* filename)
{
    std::ofstream fout(filename);
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
 * Reda a set sequence from the file and check it.
 *
 * @param filenamein: filename to read from;
 * @param filenameout: filename for output.
 */
template <int N, int D, size_t NM, size_t NMH>
void Cube_Product_Checker<N, D, NM, NMH>::read_samples_and_check(const char* filenamein, const char* filenameout)
{
    std::ifstream fin(filenamein);
    std::ofstream fout(filenameout, std::ios_base::app);
    fout << "\n ################################### \n";
    int size;
    fin >> size;
    Timewatch timer;
    double time = 0.0;
    int gv = 0;
    for (int i = 0; i < size; ++i) {
        if (N*D > 5) {
            std::cout << "working on case " << i << "\n";
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
    std::cout << "Found " << gv << " solutions\n";
    fout << "Total time: " << time << " s (" << (time/size) << " s avg) (found " << gv << " good vectors)\n";
    fout.close();
}

//=============================================================================

#endif // MMCHECKER_HPP_INCLUDED
