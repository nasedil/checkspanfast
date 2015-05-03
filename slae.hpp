/**
 * @file slae.hpp:  System of Linear Algebraic Equations utilities.
 *
 * This module contains functions for solving SLAE,
 * particularly for Cube Product Strassen-like algorithm finding.
 *
 * @author Eugene Petkevich
 * @version pre-alpha
 */

#ifndef SLAE_HPP_INCLUDED
#define SLAE_HPP_INCLUDED

#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "utils.hpp"

using namespace std;

//=============================================================================

/// All the following functions solve SLAE using different methods.

template <size_t N>
bool gauss_solve(vector<mm_bitset<N>> a, mm_bitset<N> b);

class Gauss_Presolve_Data;
template <size_t N>
void gauss_presolve(vector<mm_bitset<N>> a, Gauss_Presolve_Data& p);
template <size_t N>
bool gauss_solve(const Gauss_Presolve_Data& p, mm_bitset<N> b);

template <size_t N>
bool binary_solve(const vector<mm_bitset<N>>& a, const mm_bitset<N>& b);

template <size_t N>
bool binary_solve_recursive(const vector<mm_bitset<N>>& a, const mm_bitset<N>& b);

template <size_t N>
bool gauss_solve_randomized(vector<mm_vector_with_properties<N>> a, mm_vector_with_properties<N> b);

template <size_t N>
class Gauss_WP_Presolve_Data;
template <size_t N>
bool gauss_wp_presolve(const vector<mm_vector_with_properties<N>>& a, Gauss_WP_Presolve_Data<N>& p);
template <size_t N>
bool gauss_wp_solve(const Gauss_WP_Presolve_Data<N>& p, const mm_vector_with_properties<N>& b);

//=============================================================================

/// Comparison functions for bitvectors.

template <size_t N>
bool compare_lexical (const mm_bitset<N>& a, const mm_bitset<N>& b);

template <size_t N>
bool compare_count (const mm_bitset<N>& a, const mm_bitset<N>& b);

template <size_t N>
bool compare_combined (const mm_bitset<N>& a, const mm_bitset<N>& b);

//=============================================================================
//=============================================================================

/**
 * Sovle SLAE by Gaussian elimination.
 *
 * @param a: the coefficient matrix;
 * @param b: the resulting vector;
 *
 * @return if there is a solution.
 */
template <size_t N>
bool
gauss_solve(vector<mm_bitset<N>> a, mm_bitset<N> b)
{
    int r = 0;
    for (int i = 0; i < a.size(); ++i) { // for all columns
        // first we exchange rows if needed, so that a[i][i] = 1
        if (!a[i][r]) {
            int k = r+1;
            while ((k < N) && (!a[i][k])) {
                ++k;
            }
            if (k >= N) {
                continue;
            } else {
                bool tmp = b[k];
                b[k] = b[r];
                b[r] = tmp;
                for (int l = i; l < a.size(); ++l) {
                    tmp = a[l][k];
                    a[l][k] = a[l][r];
                    a[l][r] = tmp;
                }
            }
        }
        // second, we make all coefficients under a[i][i] equal to 0
        for (int j = r+1; j < b.size(); ++j) { // for all rows
            if (a[i][j]) {
                for (int k = i; k < a.size(); ++k) {
                    a[k][j] = (a[k][j] != a[k][r]);
                }
                b[j] = (b[j] != b[r]);
            }
        }
        ++r;
    }
    // we check all i from 7 to 15 for zero coefficients
    for (int i = r; i < b.size(); ++i) {
        if (b[i])
            return false;
    }
    return true;
}

//=============================================================================

/**
 * Data for Gaussian elimination presolve.
 */
class Gauss_Presolve_Data {
public:
    size_t n; /// Size of the matrix.
    int r;
    vector<int> rows_s;
    vector<int> rows_d;
    vector<bool> rows_o;
};

/**
 * Presolve SLAE by Gaussian elimination.
 *
 * Collects Gaussian elimination operations that are needed
 * to selve a SLAE with a specified matrix.
 *
 * @param a: the coefficient matrix;
 * @param p: varuable to store presolve data;
 */
template <size_t N>
void
gauss_presolve(vector<mm_bitset<N>> a, Gauss_Presolve_Data & p)
{
    p.n = a.size();
    p.rows_s.clear();
    p.rows_d.clear();
    p.rows_o.clear();
    int r = 0;
    for (int i = 0; i < a.size(); ++i) { // for all columns
        // first we exchange rows if needed, so that a[i][i] = 1
        if (!a[i][r]) {
            int k = r+1;
            while ((k < N) && (!a[i][k])) {
                ++k;
            }
            if (k >= N) {
                continue;
            } else {
                bool tmp;
                p.rows_s.push_back(r);
                p.rows_d.push_back(k);
                p.rows_o.push_back(true);
                for (int l = i; l < a.size(); ++l) {
                    tmp = a[l][k];
                    a[l][k] = a[l][r];
                    a[l][r] = tmp;
                }
            }
        }
        // second, we make all coefficients under a[i][i] equal to 0
        for (int j = r+1; j < N; ++j) { // for all rows
            if (a[i][j]) {
                for (int k = i; k < a.size(); ++k) {
                    a[k][j] = (a[k][j] != a[k][r]);
                }
                p.rows_s.push_back(r);
                p.rows_d.push_back(j);
                p.rows_o.push_back(false);
            }
        }
        ++r;
    }
    p.r = r;
}

/**
 * Sovle SLAE by Gaussian elimination using precomputed data.
 *
 * @param p: the presolve data;
 * @param b: the resulting vector;
 *
 * @return if there is a solution.
 */
template <size_t N>
bool
gauss_solve(const Gauss_Presolve_Data& p, mm_bitset<N> b)
{
    for (int i = 0; i < p.rows_o.size(); ++i) {
        if (p.rows_o[i]) {
            bool tmp = b[p.rows_s[i]];
            b[p.rows_s[i]] = b[p.rows_d[i]];
            b[p.rows_d[i]] = tmp;
        } else {
            b[p.rows_d[i]] = (b[p.rows_d[i]] != b[p.rows_s[i]]);
        }
    }
    for (int i = p.r; i < b.size(); ++i) {
        if (b[i])
            return false;
    }
    return true;
}

//=============================================================================

/**
 * Solve SLAE by trying all solutions.
 *
 * @param a: the coefficient matrix;
 * @param b: the resulting vector;
 *
 * @return if there is a solution.
 */
template <size_t N>
bool
binary_solve(const vector<mm_bitset<N>>& a, const mm_bitset<N>& b)
{
    uint_least64_t counter;
    uint_least64_t limit = power(2,a.size());
    for (counter = 1; counter < limit; ++counter) {
        boost::dynamic_bitset<> x(a.size(),counter);
        mm_bitset<N> r(0);
        bool is_good = true;
        // we multiply a by x
        for (int i = 0; i < b.size(); ++i) {
            for (int j = 0; j < a.size(); ++j) {
                if (x[j]) {
                    r[i] = (r[i] != a[j][i]);
                }
            }
            if (r[i] != b[i]) {
                is_good = false;
                break;
            }
        }
        if (is_good)
            return true;
    }
    return false;
}

/**
 * Solve SLAE by trying all solutions.
 *
 * @param a: the coefficient matrix;
 * @param b: the resulting vector;
 *
 * @return if there is a solution.
 */
template <size_t N>
bool
binary_solve(const vector<mm_vector_with_properties<N>>& a, const mm_vector_with_properties<N>& b)
{
    vector<mm_bitset<N>> sa;
    for (auto& vp: a) {
        sa.push_back(vp.v);
    }
    return binary_solve(sa, b.v);
}

/**
 * Solve SLAE by trying all solutions and return the result.
 *
 * The SLAE should have a solution.
 *
 * @param a: the coefficient matrix;
 * @param b: the resulting vector;
 *
 * @return solution (if SLAE has no solution, return value could be anything).
 */
template <size_t N>
boost::dynamic_bitset<>
binary_solve_result(const vector<mm_bitset<N>>& a, const mm_bitset<N>& b)
{
    uint_least64_t counter;
    uint_least64_t limit = power(2,a.size());
    for (counter = 1; counter < limit; ++counter) {
        boost::dynamic_bitset<> x(a.size(),counter);
        mm_bitset<N> r(0);
        bool is_good = true;
        // we multiply a by x
        for (int i = 0; i < b.size(); ++i) {
            for (int j = 0; j < a.size(); ++j) {
                if (x[j]) {
                    r[i] = (r[i] != a[j][i]);
                }
            }
            if (r[i] != b[i]) {
                is_good = false;
                break;
            }
        }
        if (is_good)
            return x;
    }
    return boost::dynamic_bitset<>();
}

//=============================================================================

/**
 * Subroutine for solving SLAE via backtracking.
 *
 * @param a: the coefficient matrix;
 * @param b: the resulting vector;
 * @param depth: backtracking tree current depth;
 * @param current: current value of guessed solution.
 *
 * @return if there is a solution.
 */
template <size_t N>
bool
binary_solve_explore(const vector<mm_bitset<N>>* a, const mm_bitset<N>* b, int depth, const mm_bitset<N>* current)
{
    if (depth == -1) {
        for (int k = 0; k < N; ++k) {
            if (current->test(k) != b->test(k))
                return false;
        }
        return true;
    } else {
        bool result = binary_solve_explore(a, b, depth-1, current);
        if (result)
            return true;
        mm_bitset<N>* current_add = new mm_bitset<N>();
        for (int k = 0; k < N; ++k) {
            (*current_add)[k] = current->test(k) != (*a)[depth][k];
        }
        result = binary_solve_explore(a, b, depth-1, current_add);
        delete current_add;
        return result;
    }
}

/**
 * Solve SLAE by trying solutions via backtracking.
 *
 * @param a: the coefficient matrix;
 * @param b: the resulting vector;
 *
 * @return if there is a solution.
 */
template <size_t N>
bool
binary_solve_recursive(const vector<mm_bitset<N>>& a, const mm_bitset<N>& b)
{
    uint_least64_t depth = a.size()-1;
    mm_bitset<N>* current = new mm_bitset<N>();
    current->reset();
    bool result = binary_solve_explore(&a, &b, depth, current);
    delete current;
    return result;
}

//=============================================================================

/**
 * TODO
 */
template <size_t N>
bool
gauss_solve_randomized(vector<mm_vector_with_properties<N>> a, mm_vector_with_properties<N> b)
{
    return true;
}

//=============================================================================

/**
 * Compare function to sort bitsets lexicographically.
 *
 * @param N: bitset size;
 */
template <size_t N>
bool
compare_lexical (const mm_bitset<N>& a, const mm_bitset<N>& b)
{
    for (int i = 0; i < N; ++i) {
        if (a[i] < b[i]) {
            return true;
        } else if (a[i] > b[i]) {
            return false;
        }
    }
    return true;
}

/**
 * Compare function to sort bitsets by number of 1-bits.
 *
 * @param N: bitset size;
 */
template <size_t N>
bool
compare_count (const mm_bitset<N>& a, const mm_bitset<N>& b)
{
    return (a.count() < b.count());
}

/**
 * Compare function to sort bitsets by number of bits,
 * and if it is the same, lexicographically.
 *
 * @param N: bitset size;
 */
template <size_t N>
bool
compare_combined (const mm_bitset<N>& a, const mm_bitset<N>& b)
{
    if (a.count() < b.count()) {
        return true;
    } else if (a.count() > b.count()) {
        return false;
    } else {
        for (int i = 0; i < N; ++i) {
            if (a[i] > b[i]) {
                return true;
            } else if (a[i] < b[i]) {
                return false;
            }
        }
        return true;
    }
}

//=============================================================================

/**
 * Presolve data for Gaussian elimination variant for solving SLAE.
 *
 * @param N: size of result and variable vectors.
 */
template <size_t N>
class Gauss_WP_Presolve_Data {
public:
    vector<mm_vector_with_properties<N>> am; /// Transformed Matrix.
    vector<int> imi;
};

/**
 * Collect operoationsdata for Gaussian elimination variant for solving SLAE.
 *
 * @param N: size of result and variable vectors.
 *
 * @param a: the coefficient matrix;
 * @param p: variable to store presolve data.
 */
template <size_t N>
bool
gauss_wp_presolve(const vector<mm_vector_with_properties<N>>& a, Gauss_WP_Presolve_Data<N>& p)
{
    p.am = a;
    p.imi.clear();
    for (typename vector<mm_vector_with_properties<N>>::iterator i = p.am.begin(); i != p.am.end();) { // for all columns
        int k = 0;
        while (((k < N) && (!(*i).v[k])) || (find(p.imi.begin(), p.imi.end(), k) != p.imi.end())) {
            ++k;
        }
        if (k == N) { // linearly dependent
            i = p.am.erase(i);
            return false;
        } else { // make 0's all items in a row k except for i'th column
            for (typename vector<mm_vector_with_properties<N>>::iterator j = p.am.begin(); j != i; ++j) {
                if ((*j).v[k]) {
                    (*j).v ^= (*i).v;
                    (*j).r ^= (*i).r;
                }
            }
            for (typename vector<mm_vector_with_properties<N>>::iterator j = i+1; j != p.am.end(); ++j) {
                if ((*j).v[k]) {
                    (*j).v ^= (*i).v;
                    (*j).r ^= (*i).r;
                }
            }
            p.imi.push_back(k);
            ++i;
        }
    }
    return true;
}

/**
 * Sovle SLAE by Gaussian elimination using precomputed data.
 *
 * @param N: size of result and variable vectors.
 *
 * @param p: the presolve data;
 * @param b: the resulting vector;
 *
 * @return if there is a solution.
 */
template <size_t N>
bool
gauss_wp_solve(const Gauss_WP_Presolve_Data<N>& p, const mm_vector_with_properties<N>& b)
{
    mm_vector_with_properties<N> c;
    c.r.reset();
    typename vector<mm_vector_with_properties<N>>::const_iterator j = p.am.begin();
    for (vector<int>::const_iterator i = p.imi.begin(); i != p.imi.end(); ++i) {
        if (b.v[*i]) {
            c.r ^= (*j).r;
        }
        ++j;
    }
    for (int i = 0; i < c.r.size(); ++i) {
        if (c.r[i] != b.r[i])
            return false;
    }
    c.v.reset();
    j = p.am.begin();
    for (vector<int>::const_iterator i = p.imi.begin(); i != p.imi.end(); ++i) {
        if (b.v[*i]) {
            c.v ^= (*j).v;
        }
        ++j;
    }
    for (int i = 0; i < c.v.size(); ++i) {
        if (c.v[i] != b.v[i])
            return false;
    }
    return true;
}

//=============================================================================

#endif // SLAE_HPP_INCLUDED
