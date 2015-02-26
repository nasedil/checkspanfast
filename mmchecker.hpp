#ifndef MMCHECKER_HPP_INCLUDED
#define MMCHECKER_HPP_INCLUDED

#include <algorithm>
#include <iostream>
#include <fstream>
#include <set>
#include <vector>

#include "slae.hpp"
#include "utils.hpp"

//=======================================================

/**
 * Main class, calculates all nesessary information
 * for given template parameters:
 *
 * N: size of the cube;
 * D: dimension of the cube;
 * NM: number of bits in multiplication vectors ((N^D)^D);
 * NMH: number of elements in the array (N^D).
 */
template <int N, int D, size_t NM, size_t NMH>
class Matrix_Multiplication_Checker
{
public:
	typedef mm_bitset<NM> Multiplication_Vector; // multiplication vector type
	typedef mm_bitset<NMH> Multiplication_Part_Vector; // element coefficient type
    int length; // size of cube (N)
    int dimension; // dimension of cube (D)
    int element_count; // number of elements in cube (N^D)
    mm_vector_with_properties<NM> * r_vectors; // result vectors
    mm_vector_with_properties<NM> * m_vectors; // multiplication vectors
    int m_count; // number of multiplication vectors = 2^((N^D)*D)
    int m_length; // number of non-zero linear combinations = 2^((N^D-1)*D)
    int f_count; // number of vectors to choose for a basis = N^(D+1)-1
    std::vector<Multiplication_Vector> good_vectors; // set of vectors that are in the current span
    std::set<int> good_vectors_indexes; // set of indexes for these vectors
    std::vector<Multiplication_Vector> n_vectors; // set of current chosen vectors
    std::set<int> n_vectors_indexes; // set of their indexes
    Timewatch tw; // timer that is used for getting calculation time

    //=============--- constructors and destructors
    Matrix_Multiplication_Checker();
    ~Matrix_Multiplication_Checker();

    //=============--- index operations
    // return linear index in a bit vector by its indexes in matrices
    int get_vector_index(int ai, int aj, int bi, int bj); // for 2-dimensianal case
    int get_vector_index(int ai, int aj, int ak, int bi, int bj, int bk, int ci, int cj, int ck); // for 3-d case
    // return linear index in a bit vector by combined indexes in matrices
    int get_vector_index(int a, int b); // for 2-d case
    int get_vector_index(int a, int b, int c); // for 3-d case
	// return indices from index
    void decode_indices_from_index(int index, int & ai, int & aj, int & bi, int & bj); // for 2-d case
    void decode_indices_from_index(int index, int & ai, int & aj, int & ak, int & bi, int & bj, int & bk, int & ci, int & cj, int & ck); // for 3-d case
    // return linear index of an element in a matrix
    int get_element_index(int i, int j); // for 2-d case
    int get_element_index(int i, int j, int k); // for 3-d case
    // return index of vector in m
    int get_m_index(int i, int j); // for 2-d case
    int get_m_index(int i, int j, int k); // for 3-d case

    //=============--- Initial calculations
    // write result vectors to array
    void calculate_r_vectors();
    // write multiplication vectors to array
    void calculate_m_vectors();

    //=============--- checking routines
    void clear_sets(); // clear current sets of vectors
    void add_vector_to_set(int cc); // add vector to current sets
    bool check_vectors_for_goodness(); // check current set of vectors

    //=============--- searching for solution
    bool check_for_good_vectors(); // check all solution space
    bool check_for_good_vectors_randomized(); // do random walk

    //=============--- utilities
    void output_vector(Multiplication_Vector v); // output vector to screen
    void output_vector_text(Multiplication_Vector v); // output vector to screen in letters
    void save_random_samples(int size, const char * filename); // save random sets to a file (for testing later)
    void read_samples_and_check(const char * filename, const char * filenameout); // check sets from a file
    bool check_vectors_for_goodness_a1(); // just normal gauss everywhere
};

//=======================================================
//=======================================================

template <int N, int D, size_t NM, size_t NMH>
Matrix_Multiplication_Checker<N, D, NM, NMH>::Matrix_Multiplication_Checker()
{
	std::cout << "MMChecker created.\n";
    this->length = N;
    this->dimension = D;
    this->element_count = power(this->length, this->dimension);
    this->f_count = power(length,dimension+1)-1;
    mm_vector_with_properties_options o;
    mm_vector_with_properties<NM>::make_options(o);
    // calculate r_vectors
    tw.watch();
    r_vectors = new mm_vector_with_properties<NM>[element_count];
    calculate_r_vectors();
    for (int i = 0; i < element_count; ++i)
	{
		r_vectors[i].calculate_properties(o);
	}
    std::cout << "[" << tw.watch() << " s] Result vectors calculated.\n";
    // calculate m_vectors
    m_length = power(2,element_count)-1;
    m_count = power(m_length,dimension);
    m_vectors = new mm_vector_with_properties<NM>[m_count];
    calculate_m_vectors();
    for (int i = 0; i < m_count; ++i)
	{
		m_vectors[i].calculate_properties(o);
	}
    std::cout << "[" << tw.watch() << " s] Multiplication vectors calculated.\n";
}

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
Matrix_Multiplication_Checker<N, D, NM, NMH>::~Matrix_Multiplication_Checker()
{
    delete [] r_vectors;
    delete [] m_vectors;
}

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
inline int Matrix_Multiplication_Checker<N, D, NM, NMH>::get_vector_index(int ai, int aj, int bi, int bj)
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

template <int N, int D, size_t NM, size_t NMH>
inline int Matrix_Multiplication_Checker<N, D, NM, NMH>::get_vector_index(int ai, int aj, int ak, int bi, int bj, int bk, int ci, int cj, int ck)
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

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
inline int Matrix_Multiplication_Checker<N, D, NM, NMH>::get_vector_index(int a, int b)
{
    int result = a;
    result *= element_count;
    result += b;
    return result;
}

template <int N, int D, size_t NM, size_t NMH>
inline int Matrix_Multiplication_Checker<N, D, NM, NMH>::get_vector_index(int a, int b, int c)
{
	int result = a;
    result *= element_count;
    result += b;
    result *= element_count;
    result += c;
    return result;
}

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
inline void Matrix_Multiplication_Checker<N, D, NM, NMH>::decode_indices_from_index(int index, int & ai, int & aj, int & bi, int & bj)
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

template <int N, int D, size_t NM, size_t NMH>
inline void Matrix_Multiplication_Checker<N, D, NM, NMH>::decode_indices_from_index(int index, int & ai, int & aj, int & ak, int & bi, int & bj, int & bk, int & ci, int & cj, int & ck)
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

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
inline int Matrix_Multiplication_Checker<N, D, NM, NMH>::get_element_index(int i, int j)
{
    return ((i*length) + j);
}

template <int N, int D, size_t NM, size_t NMH>
inline int Matrix_Multiplication_Checker<N, D, NM, NMH>::get_element_index(int i, int j, int k)
{
	return (((i*length) + j)*length + k);
}

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
inline int Matrix_Multiplication_Checker<N, D, NM, NMH>::get_m_index(int i, int j)
{
    return ((i*m_length) + j);
}

template <int N, int D, size_t NM, size_t NMH>
inline int Matrix_Multiplication_Checker<N, D, NM, NMH>::get_m_index(int i, int j, int k)
{
	return (((i*m_length) + j)*m_length + k);
}

//=======================================================

template <>
void Matrix_Multiplication_Checker<2, 2, 16, 4>::calculate_r_vectors()
{
    for (int i = 0; i < length; ++i)
        for (int j = 0; j < length; ++j)
        {
            int index = get_element_index(i, j);
            r_vectors[index].v.reset();
            for (int l = 0; l < length; ++l)
            {
                r_vectors[index].v[get_vector_index(i, l, l, j)] = 1;
            }
        }
}

template <>
void Matrix_Multiplication_Checker<2, 3, 512, 8>::calculate_r_vectors()
{
    for (int i = 0; i < length; ++i)
        for (int j = 0; j < length; ++j)
			for (int k = 0; k < length; ++k)
			{
				int index = get_element_index(i, j, k);
				r_vectors[index].v.reset();
				for (int l = 0; l < length; ++l)
				{
					r_vectors[index].v[get_vector_index(i, j, l, i, l, k, l, j, k)] = 1;
				}
			}
}

//=======================================================

template <>
void Matrix_Multiplication_Checker<2, 2, 16, 4>::calculate_m_vectors()
{
    for (int i = 1; i < power(2,element_count); ++i)
        for (int j = 1; j < power(2,element_count); ++j)
		{
			Multiplication_Part_Vector av(i);
			Multiplication_Part_Vector bv(j);
			int index = get_m_index(i-1, j-1);
			m_vectors[index].v.reset();
			for (int k = 0; k < element_count; ++k)
				for (int l = 0; l < element_count; ++l)
				{
					if (av[k])
						if (bv[l])
							m_vectors[index].v.set(get_vector_index(k, l));
				}
		}
}

template <>
void Matrix_Multiplication_Checker<2, 3, 512, 8>::calculate_m_vectors()
{
    for (int i = 1; i < power(2,element_count); ++i)
        for (int j = 1; j < power(2,element_count); ++j)
			for (int k = 1; k < power(2,element_count); ++k)
			{
				//cout << i << " " << j << " " << k << "\n";
				Multiplication_Part_Vector av(i);
				Multiplication_Part_Vector bv(j);
				Multiplication_Part_Vector cv(k);
				int index = get_m_index(i-1, j-1, k-1);
				m_vectors[index].v.reset();
				for (int l = 0; l < element_count; ++l)
					for (int o = 0; o < element_count; ++o)
						for (int p = 0; p < element_count; ++p)
						{
							if (av[l])
								if (bv[o])
									if (cv[p])
										m_vectors[index].v.set(get_vector_index(l, o, p));
						}
			}
}

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
bool Matrix_Multiplication_Checker<N, D, NM, NMH>::check_vectors_for_goodness()
{
#ifdef output1
	std::cout << "Checking of vectors { ";
	for (int i: n_vectors_indexes)
		std::cout << i << " ";
	std::cout << "} started.\n";
#endif
	std::vector<mm_vector_with_properties<NM>> nvwp;
	std::vector<mm_vector_with_properties<NM>> gvwp;
	for (std::set<int>::iterator i = n_vectors_indexes.begin(); i != n_vectors_indexes.end(); ++i)
		nvwp.push_back(m_vectors[*i]);
	for (std::set<int>::iterator i = good_vectors_indexes.begin(); i != good_vectors_indexes.end(); ++i)
		gvwp.push_back(m_vectors[*i]);
	for (int i = 0; i < element_count; ++i)
	{
		n_vectors.push_back(r_vectors[i].v);
		nvwp.push_back(r_vectors[i]);
	}
	Vectors_Presolve_Data<NM> v;
	for (int i = 0; i < n_vectors.size(); ++i)
		v.add_vector(n_vectors[i]);
	Gauss_WP_Presolve_Data<NM> pwp;
	gauss_wp_presolve(nvwp, pwp);
	Gauss_WP_Presolve_Data<NM> pwpg;
	for (int i = 0; i < m_count; ++i)
	{
		if (!v.check(m_vectors[i].v))
			continue;
		if (gauss_wp_solve(pwp, m_vectors[i]))
		{
			gauss_wp_presolve(gvwp, pwpg);
			if (!gauss_wp_solve(pwpg, m_vectors[i]))
			{
				good_vectors.push_back(m_vectors[i].v);
				good_vectors_indexes.insert(i);
				gvwp.push_back(m_vectors[i]);
			}
			if (good_vectors.size() >= f_count)
				break;
		}
	}
	if (good_vectors.size() >= f_count)
	{
		std::cout << "Good vectors found: ";
		for (std::set<int>::iterator cc = good_vectors_indexes.begin(); cc != good_vectors_indexes.end(); ++cc)
			std::cout << *cc << " ";
		std::cout << "\n";
#ifdef output2
		std::ofstream mvfile("good.txt");
		mvfile << "Found good vectors!!!\n";
		for (std::set<int>::iterator cc = good_vectors_indexes.begin(); cc != good_vectors_indexes.end(); ++cc)
			mvfile << *cc << " ";
		mvfile << "\n";
#endif
		return true;
	}
	return false;
}

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
void Matrix_Multiplication_Checker<N, D, NM, NMH>::clear_sets()
{
	n_vectors.clear();
	n_vectors_indexes.clear();
	good_vectors.clear();
	good_vectors_indexes.clear();
}

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
void Matrix_Multiplication_Checker<N, D, NM, NMH>::add_vector_to_set(int cc)
{
	n_vectors.push_back(m_vectors[cc].v);
	n_vectors_indexes.insert(cc);
	//good_vectors.push_back(m_vectors[cc].v);
	//good_vectors_indexes.insert(cc);
}

//=======================================================

template <>
bool Matrix_Multiplication_Checker<2, 2, 16, 4>::check_for_good_vectors()
{
	std::set<std::set<int>> results;
    for (int c1 = 0; c1 < m_count-2; ++c1)
	{
		std::cout << c1 << "\n";
		for (int c2 = c1+1; c2 < m_count-1; ++c2)
		{
			for (int c3 = c2+1; c3 < m_count; ++c3)
			{
				// 9  37  62  109  128  146  165
				clear_sets();
				add_vector_to_set(c1);
				add_vector_to_set(c2);
				add_vector_to_set(c3);
				if (check_vectors_for_goodness())
					results.insert(good_vectors_indexes);
			}
		}
	}
	std::ofstream fout("results2.txt");
	fout << "Results: " << results.size() << " found solutions.\n";
	for (std::set<std::set<int>>::iterator i = results.begin(); i != results.end(); ++i)
	{
		fout << "[ ";
		for (std::set<int>::iterator j = (*i).begin(); j != (*i).end(); ++j)
		{
			fout << (*j) << " ";
		}
		fout << "]\n";
	}
	fout.close();
	return false;
}

template <>
bool Matrix_Multiplication_Checker<2, 3, 512, 8>::check_for_good_vectors()
{
    for (int c1 = 0; c1 < m_count-6; ++c1)
		for (int c2 = c1+1; c2 < m_count-5; ++c2)
			for (int c3 = c2+1; c3 < m_count-4; ++c3)
				for (int c4 = c3+1; c4 < m_count-3; ++c4)
					for (int c5 = c4+1; c5 < m_count-2; ++c5)
						for (int c6 = c5+1; c6 < m_count-1; ++c6)
							for (int c7 = c6+1; c7 < m_count; ++c7)
							{
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

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
bool Matrix_Multiplication_Checker<N, D, NM, NMH>::check_for_good_vectors_randomized()
{
	std::ofstream mvfile("mv.txt", std::ios_base::app);
	// TODO read and write file
	Random rnd(0, m_count-1);
	while (true)
	{
		clear_sets();
		for (int i = 0; i < (f_count-element_count); ++i)
		{
			int cc = rnd.next();
			while (n_vectors_indexes.count(cc) > 0)
			{
				cc = rnd.next();
			}
			add_vector_to_set(cc);
		}
		if (check_vectors_for_goodness())
			return true;
		else
		{
			// TODO write to file (or memory first)
		}
	}
	return false;
}

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
void Matrix_Multiplication_Checker<N, D, NM, NMH>::output_vector(Multiplication_Vector v)
{
	std::cout << v;
}

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
void Matrix_Multiplication_Checker<N, D, NM, NMH>::output_vector_text(Multiplication_Vector v)
{
    for (int i = 0; i < v.size(); ++i)
    {
        if (v[i])
        {
            int ai, aj, bi, bj;
            decode_indices_from_index(i, ai, aj, bi, bj);
            std::cout << "A" << ai+1 << aj+1 << "B" << bi+1 << bj+1 << " ";
        }
    }
}

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
void Matrix_Multiplication_Checker<N, D, NM, NMH>::save_random_samples(int size, const char * filename)
{
	std::ofstream fout(filename);
	Random rnd(0, m_count-1);
	fout << size << "\n";
	for (int i = 0; i < size; ++i)
	{
		clear_sets();
		for (int i = 0; i < (f_count-element_count); ++i)
		{
			int cc = rnd.next();
			while (n_vectors_indexes.count(cc) > 0)
			{
				cc = rnd.next();
			}
			add_vector_to_set(cc);
			fout << cc << " ";
		}
		fout << "\n";
	}
	fout.close();
}

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
void Matrix_Multiplication_Checker<N, D, NM, NMH>::read_samples_and_check(const char * filenamein, const char * filenameout)
{
	std::ifstream fin(filenamein);
	std::ofstream fout(filenameout, std::ios_base::app);
	fout << "\n ################################### \n";
	int size;
	fin >> size;
	Timewatch timer;
	double time = 0.0;
	int gv = 0;
	for (int i = 0; i < size; ++i)
	{
		if (N*D > 5)
			std::cout << "working on case " << i << "\n";
		clear_sets();
		for (int i = 0; i < (f_count-element_count); ++i)
		{
			int cc;
			fin >> cc;
			add_vector_to_set(cc);
		}
		timer.watch();
		if (check_vectors_for_goodness())
			++gv;
		double curtime = timer.watch();
		time += curtime;
		if (N*D > 5)
			fout << "\t" << i << ": " << curtime << " s\n";
	}
	std::cout << "Found " << gv << " solutions\n";
	fout << "Total time: " << time << " s (" << (time/size) << " s avg) (found " << gv << " good vectors)\n";
	fout.close();
}

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
bool Matrix_Multiplication_Checker<N, D, NM, NMH>::check_vectors_for_goodness_a1()
{
	for (int i = 0; i < element_count; ++i)
	{
		n_vectors.push_back(r_vectors[i].v);
	}
	std::sort(n_vectors.begin(), n_vectors.end(), compare_combined<NM>);
	for (int i = 0; i < m_count; ++i)
	{
		if (gauss_solve(n_vectors, m_vectors[i].v))
		{
			if (!gauss_solve(good_vectors, m_vectors[i].v))
			{
				good_vectors.push_back(m_vectors[i].v);
				good_vectors_indexes.insert(i);
				std::sort(good_vectors.begin(), good_vectors.end(), compare_combined<NM>);
			}
			if (good_vectors.size() >= f_count)
				break;
		}
	}
	if (good_vectors.size() >= f_count)
		return true;
	return false;
}

//=======================================================

#endif // MMCHECKER_HPP_INCLUDED
