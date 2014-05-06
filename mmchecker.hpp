#ifndef MMCHECKER_HPP_INCLUDED
#define MMCHECKER_HPP_INCLUDED

#include <ctime>
#include <algorithm>
#include <fstream>
#include <set>
#include "slae.hpp"
#include "utils.hpp"

template <int N, int D, size_t NM, size_t NMH>
class Matrix_Multiplication_Checker
{
public:
	typedef mm_bitset<NM> Multiplication_Vector;
	typedef mm_bitset<NMH> Multiplication_Part_Vector;
    int length;
    int dimension;
    int element_count;
    mm_vector_with_properties<NM> * r_vectors; // result vectors
    mm_vector_with_properties<NM> * m_vectors; // multiplication vectors
    int m_count; // number of multiplication vectors
    int m_length; // number of non-zero linear combinations
    int f_count; // number of vectors to choose for a basis
    vector<Multiplication_Vector> good_vectors;
    set<int> good_vectors_indexes;
    vector<Multiplication_Vector> n_vectors;
    set<int> n_vectors_indexes;

    Timewatch tw;

    Matrix_Multiplication_Checker();
    ~Matrix_Multiplication_Checker();

    // returns linear index in a bit vector by its indexes in matrices
    int get_vector_index(int ai, int aj, int bi, int bj);
    int get_vector_index(int ai, int aj, int ak, int bi, int bj, int bk, int ci, int cj, int ck);
    // returns linear index in a bit vector by combined indexes in matrices
    int get_vector_index(int a, int b);
    int get_vector_index(int a, int b, int c);
	// returns indices from index
    void decode_indices_from_index(int index, int & ai, int & aj, int & bi, int & bj);
    void decode_indices_from_index(int index, int & ai, int & aj, int & ak, int & bi, int & bj, int & bk, int & ci, int & cj, int & ck);
    // returns linear index of an element in a matrix
    int get_element_index(int i, int j);
    int get_element_index(int i, int j, int k);
    // returns index of vector in m
    int get_m_index(int i, int j);
    int get_m_index(int i, int j, int k);
    // write result vectors to array
    void calculate_r_vectors();
    // write different multiplication vectors to array
    void calculate_m_vectors();
    // outputs vector to screen
    void output_vector(Multiplication_Vector v);
    // outputs vector to screen in letters
    void output_vector_text(Multiplication_Vector v);
    // finds if there are good vectors (what we want to find)
    void clear_sets();
    void add_vector_to_set(int cc);
    bool check_for_good_vectors();
    bool check_for_good_vectors_randomized();
    bool check_vectors_for_goodness();
    void cvfg_precalculate();
    bool cvfg_slae();
};

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
Matrix_Multiplication_Checker<N, D, NM, NMH>::Matrix_Multiplication_Checker()
{
	cout << "MMChecker created.\n";
    this->length = N;
    this->dimension = D;
    this->element_count = power(this->length, this->dimension);
    this->f_count = power(length,dimension+1)-1;
    mm_vector_with_properties_options o;
    mm_vector_with_properties<NM>::make_options(o);
    /*
    for (int i = 0; i < o.rnd_sets.size(); ++i)
	{
		cout << "Set #" << i << ": ";
		for (set<int>::iterator j = o.rnd_sets[i].begin(); j != o.rnd_sets[i].end(); ++j)
			cout << (*j) << " ";
		cout << "\n";
	}
	*/
    // calculate r_vectors
    tw.watch();
    r_vectors = new mm_vector_with_properties<NM>[element_count];
    calculate_r_vectors();
    for (int i = 0; i < element_count; ++i)
	{
		r_vectors[i].calculate_properties(o);
		//cout << "Result vector #" << i << " : " << r_vectors[i].v << " " << r_vectors[i].r << "\n";
	}
    cout << "[" << tw.watch() << " s] Result vectors calculated.\n";
    // calculate m_vectors
    m_length = power(2,element_count)-1;
    m_count = power(m_length,dimension);
    m_vectors = new mm_vector_with_properties<NM>[m_count];
    calculate_m_vectors();
    for (int i = 0; i < m_count; ++i)
	{
		m_vectors[i].calculate_properties(o);
		//cout << "Mult. vector #" << i << " : " << m_vectors[i].v << " " << m_vectors[i].r << "\n";
	}
    cout << "[" << tw.watch() << " s] Multiplication vectors calculated.\n";
    // additional calculations
    //sort(m_vectors, m_vectors+m_count, compare_combined<NM>);
    //cout << "[" << tw.watch() << " s] Multiplication vectors sorted.\n";

	//ofstream mvfile("mv.txt");
	//for (int i = 0; i < m_count; ++i)
	//	mvfile << m_vectors[i] << "\n";
	//mvfile.close();
    //cout << "[" << tw.watch() << " s] Multiplication vectors stored.\n";
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
				//cout << m_vectors[index] << "\n";
			}
}

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
bool Matrix_Multiplication_Checker<N, D, NM, NMH>::check_vectors_for_goodness()
{
	tw.watch();
	//cout << ".";
#ifdef output1
	cout << "Checking of vectors { ";
	for (int i: n_vectors_indexes)
		cout << i << " ";
	cout << "} started.\n";
#endif
	//tw.watch();
	//sort(good_vectors.begin(), good_vectors.end(), compare_combined<NM>);
	vector<mm_vector_with_properties<NM>> nvwp;
	vector<mm_vector_with_properties<NM>> gvwp;
	for (set<int>::iterator i = n_vectors_indexes.begin(); i != n_vectors_indexes.end(); ++i)
		nvwp.push_back(m_vectors[*i]);
	for (set<int>::iterator i = good_vectors_indexes.begin(); i != good_vectors_indexes.end(); ++i)
		gvwp.push_back(m_vectors[*i]);
	for (int i = 0; i < element_count; ++i)
	{
		n_vectors.push_back(r_vectors[i].v);
		nvwp.push_back(r_vectors[i]);
	}
	//sort(n_vectors.begin(), n_vectors.end(), compare_combined<NM>);
	//cout << "[" << tw.watch() << " s] Vectors sorted.\n";
	// gauss presolve
	//Gauss_Presolve_Data p;
	//gauss_presolve(n_vectors,p);
	Vectors_Presolve_Data<NM> v;
	for (int i = 0; i < n_vectors.size(); ++i)
		v.add_vector(n_vectors[i]);
	// randomized gauss presolve
	Gauss_WP_Presolve_Data<NM> pwp;
	gauss_wp_presolve(nvwp, pwp);
	Gauss_WP_Presolve_Data<NM> pwpg;
	/*
	cout << "Indexes : ";
	for (int i = 0; i < pwp.imi.size(); ++i)
		cout << (pwp.imi[i]) << " ";
	cout << "\n";
	for (int i = 0; i < pwp.am.size(); ++i)
		cout << "N vector #" << i << " : " << pwp.am[i].v << " " << pwp.am[i].r << "\n";
	*/
	//cout << "[" << tw.watch() << " s] Gauss presolve done(" << p.rows_o.size() << " ops).\n";
#ifdef output1
	//cout << "Gauss presolve done(" << p.rows_o.size() << " ops).\n";
#endif
	int n_index = 0;
	for (int i = n_index; i < m_count; ++i)
	{
		// TODO: make a tree query here
		/*
		if (n_vectors_indexes.count(i) > 0)
			continue;
			*/
		//tw.watch();
		if (!v.check(m_vectors[i].v))
			continue;
		//if (gauss_solve(p, m_vectors[i].v))
		if (gauss_wp_solve(pwp, m_vectors[i]))
		{
			//cout << "[" << tw.watch() << " s] Potential vector (" << (i+1) << " from " << m_count << ") found.\n";
			//tw.watch();
			gauss_wp_presolve(gvwp, pwpg);
			//if (!gauss_solve(good_vectors,m_vectors[i].v))
			if (!gauss_wp_solve(pwpg, m_vectors[i]))
			{
				good_vectors.push_back(m_vectors[i].v);
				good_vectors_indexes.insert(i);
				gvwp.push_back(m_vectors[i]);
				//sort(good_vectors.begin(), good_vectors.end(), compare_combined<NM>);
			}
			if (good_vectors.size() >= f_count)
				break;
			//cout << "[" << tw.watch() << " s] Gauss elimination check done.\n";
		}
		else
		{
			//cout << "[" << tw.watch() << " s] Vector (" << (i+1) << " from " << m_count << ") discarded.\n";
		}
	}
	if (good_vectors.size() >= f_count)
	{
		cout << "Good vectors found: ";
		for (set<int>::iterator cc = good_vectors_indexes.begin(); cc != good_vectors_indexes.end(); ++cc)
			cout << *cc << " ";
		cout << "\n";

#ifdef output2
		ofstream mvfile("good.txt");
		mvfile << "Found good vectors!!!\n";
		for (set<int>::iterator cc = good_vectors_indexes.begin(); cc != good_vectors_indexes.end(); ++cc)
			mvfile << *cc << " ";
		mvfile << "\n";
#endif
		return true;
	}
#ifdef output1
	cout << "[" << tw.watch() << " s] Checkng finished.\n";
#endif
	return false;
}

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
void Matrix_Multiplication_Checker<N, D, NM, NMH>::cvfg_precalculate()
{
	return;
}

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
bool Matrix_Multiplication_Checker<N, D, NM, NMH>::cvfg_slae()
{
	return true;
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
	set<set<int>> results;
    for (int c1 = 0; c1 < m_count-2; ++c1)
	{
		cout << c1 << "\n";
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
	ofstream fout("results2.txt");
	fout << "Results: " << results.size() << " found solutions.\n";
	for (set<set<int>>::iterator i = results.begin(); i != results.end(); ++i)
	{
		fout << "[ ";
		for (set<int>::iterator j = (*i).begin(); j != (*i).end(); ++j)
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
	ofstream mvfile("mv.txt", std::ios_base::app);
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
	cout << v;
	/*
    for (int i = 0; i < v.size(); ++i)
    {
        cout << v[i];
    }
    */
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
            cout << "A" << ai+1 << aj+1 << "B" << bi+1 << bj+1 << " ";
        }
    }
}

//=======================================================

#endif // MMCHECKER_HPP_INCLUDED
