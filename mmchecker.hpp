#ifndef MMCHECKER_HPP_INCLUDED
#define MMCHECKER_HPP_INCLUDED

#include "slae.hpp"

template <int N, int D, size_t NM, size_t NMH>
class Matrix_Multiplication_Checker
{
public:
	typedef bitset<NM> Multiplication_Vector;
	typedef bitset<NMH> Multiplication_Part_Vector;
    int length;
    int dimension;
    int element_count;
    Multiplication_Vector * r_vectors; // result vectors
    Multiplication_Vector * m_vectors; // multiplication vectors
    int m_count; // number of multiplication vectors
    int m_length; // number of non-zero linear combinations
    int f_count; // number of vectors to choose for a basis
    vector<Multiplication_Vector> good_vectors;
    vector<int> good_vectors_indexes;

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
    bool check_for_good_vectors();
};

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
Matrix_Multiplication_Checker<N, D, NM, NMH>::Matrix_Multiplication_Checker()
{
    this->length = N;
    this->dimension = D;
    this->element_count = power(this->length, this->dimension);
    this->f_count = power(length,dimension+1)-1;
    // calculate r_vectors
    r_vectors = new Multiplication_Vector[element_count];
    calculate_r_vectors();
    // calculate m_vectors
    m_length = power(2,element_count)-1;
    m_count = power(m_length,dimension);
    m_vectors = new Multiplication_Vector[m_count];
    calculate_m_vectors();
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
            r_vectors[index].reset();
            for (int l = 0; l < length; ++l)
            {
                r_vectors[index][get_vector_index(i, l, l, j)] = 1;
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
				r_vectors[index].reset();
				for (int l = 0; l < length; ++l)
				{
					r_vectors[index][get_vector_index(i, j, l, i, l, k, l, j, k)] = 1;
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
			m_vectors[index].reset();
			for (int k = 0; k < element_count; ++k)
				for (int l = 0; l < element_count; ++l)
				{
					if (av[k])
						if (bv[l])
							m_vectors[index].set(get_vector_index(k, l));
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
				Multiplication_Part_Vector av(i);
				Multiplication_Part_Vector bv(j);
				Multiplication_Part_Vector cv(k);
				int index = get_m_index(i-1, j-1, k-1);
				m_vectors[index].reset();
				for (int l = 0; l < element_count; ++l)
					for (int o = 0; o < element_count; ++o)
						for (int p = 0; p < element_count; ++p)
						{
							if (av[l])
								if (bv[o])
									if (cv[p])
										m_vectors[index].set(get_vector_index(l, o, p));
						}
			}
}

//=======================================================

template <int N, int D, size_t NM, size_t NMH>
bool Matrix_Multiplication_Checker<N, D, NM, NMH>::check_for_good_vectors()
{
	int gvcounter = 0;
	clock_t timestamp = clock();
	int_least64_t go_counter = 0;
	int slau_counter = 0;
	int slau_counter2 = 0;
	int slau_counter3 = 0;
    sort(m_vectors, m_vectors+m_count, compare_combined<16>);
    for (int c1 = 0; c1 < m_count-2; ++c1)
	{
		double elapsed_secs = double(clock() - timestamp) / CLOCKS_PER_SEC;
		timestamp = clock();
		cout << "time: " << elapsed_secs << "\n";
		cout << "foundv: " << gvcounter << "\n";
		gvcounter = 0;
		//cout << "\nC1 = " << c1 << "\n";
		for (int c2 = c1+1; c2 < m_count-1; ++c2)
		{
			//cout << "c=" << c2;
			for (int c3 = c2+1; c3 < m_count; ++c3)
			{
				//cout << ".";
				//cout << "\t" << c1 << "\t" << c2 << "\t" << c3 << "\n";
				good_vectors.clear();
				good_vectors_indexes.clear();
				vector<Multiplication_Vector> n_vectors;
				// 9  37  62  109  128  146  165
				n_vectors.push_back(m_vectors[c1]);
				n_vectors.push_back(m_vectors[c2]);
				n_vectors.push_back(m_vectors[c3]);
				good_vectors.push_back(m_vectors[c1]);
				good_vectors.push_back(m_vectors[c2]);
				good_vectors.push_back(m_vectors[c3]);
				good_vectors_indexes.push_back(c1);
				good_vectors_indexes.push_back(c2);
				good_vectors_indexes.push_back(c3);
				sort(good_vectors.begin(), good_vectors.end(), compare_combined<16>);
				for (int i = 0; i < element_count; ++i)
					n_vectors.push_back(r_vectors[i]);
				int n_index = 0;
				sort(n_vectors.begin(), n_vectors.end(), compare_combined<16>);
				//reverse(n_vectors.begin(), n_vectors.end());
				Gauss_Presolve_Data p;
				gauss_presolve(n_vectors,p);
				// statistics
				slau_counter++;
				go_counter += p.rows_o.size();
				// end statistics
				for (int i = n_index; i < m_count; ++i)
				{
					if ((i == c1) || (i == c2) || (i == c3))
						continue;
					++slau_counter2;
					if (gauss_solve(p, m_vectors[i]))
					{
						//if (!binary_solve_recursive(n_vectors, m_vectors[i]))
						//	cout << "Error in solvers!\n";
						if (good_vectors.size() == 0)
						{
							good_vectors.push_back(m_vectors[i]);
							good_vectors_indexes.push_back(i);
						}
						else
						{
							++slau_counter3;
							sort(good_vectors.begin(), good_vectors.end(), compare_combined<16>);
							if (!gauss_solve(good_vectors,m_vectors[i]))
							{
								good_vectors.push_back(m_vectors[i]);
								good_vectors_indexes.push_back(i);
							}
						}
						if (good_vectors.size() >= f_count)
							break;
						//cout << "number " << i << "\n";
						//output_vector_text(m_vectors[i]);
						//cout << "\n";
					}
				}
				if (good_vectors.size() >= f_count)
				{
					cout << "Found good vectors (" << c1 << " " << c2 << " " << c3 << "): ";
					for (int cc = 0; cc < good_vectors.size(); ++cc)
						cout << good_vectors_indexes[cc] << " ";
					cout << "\n";
					gvcounter++;
				}
			}
		}
	}
	long double avg_go = go_counter / (double) slau_counter;
	cout << "\nAvg operations in presolve: " << avg_go << "\n";
	cout << "\nSlae solved 1: " << slau_counter2 << "\n";
	cout << "\nSlae solved 2: " << slau_counter3 << "\n";
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
