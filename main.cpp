#include <iostream>
#include <bitset>
#include <vector>
#include <fstream>
#include <ctime>
#include <algorithm>
#include "types.hpp"
#include "slae.hpp"
using namespace std;

// for positive b
int power(int a, int b);

class Matrix_Multiplication_Checker
{
public:
    int length;
    int element_count;
    Bitvector_16 * r_vectors;
    Bitvector_16 * m_vectors;
    int m_count;
    int m_length;
    int f_count;
    vector<Bitvector_16> good_vectors;
    vector<int> good_vectors_indexes;

    Matrix_Multiplication_Checker(int length);
    ~Matrix_Multiplication_Checker();

    // returns linear index in a bit vector by 4 indexes in matrices
    int get_vector_index(int ai, int aj, int bi, int bj);
    // returns linear index in a bit vector by 2 indexes in matrices
    int get_vector_index(int a, int b);
    // returns 4 indices from index
    void decode_indices_from_index(int index, int & ai, int & aj, int & bi, int & bj);
    // returns linear index of an elementt in a matrix
    int get_element_index(int i, int j);
    // returns index of vertar in m
    int get_m_index(int i, int j);
    // write result vectors to array
    void calculate_r_vectors();
    // write different multiplication vectors to array
    void calculate_m_vectors();
    // outputs vector to screen
    void output_vector(Bitvector_16 v);
    // outputs vector to screen in letters
    void output_vector_text(Bitvector_16 v);
    // finds if there are good vectors (what we want to find)
    bool check_for_good_vectors();
};

int main()
{
    // generate all possiple multiplications (> 200 vectors)

    Matrix_Multiplication_Checker * checker = new Matrix_Multiplication_Checker(2);

    /*
    for (int i = 0; i < checker->length; ++i)
        for (int j = 0; j < checker->length; ++j)
        {
            checker->output_vector_text(checker->r_vectors[checker->get_element_index(i, j)]);
            cout << "\n";
        }

    for (int i = 0; i < checker->m_count; ++i)
    {
        cout << i << " : ";
        checker->output_vector_text(checker->m_vectors[i]);
        cout << "\n";
    }
    */

    checker->check_for_good_vectors();

    delete checker;

    return 0;
}

Matrix_Multiplication_Checker::Matrix_Multiplication_Checker(int length)
{
    this->length = length;
    this->element_count = this->length * this->length;
    // calculate r_vectors
    r_vectors = new Bitvector_16[length*length];
    calculate_r_vectors();
    m_length = power(2,element_count)-1;
    m_count = m_length*m_length;
    m_vectors = new Bitvector_16[m_count];
    calculate_m_vectors();
    f_count = length*length*length-1;
}

Matrix_Multiplication_Checker::~Matrix_Multiplication_Checker()
{
    delete [] r_vectors;
    delete [] m_vectors;
}

inline int Matrix_Multiplication_Checker::get_vector_index(int ai, int aj, int bi, int bj)
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

inline void Matrix_Multiplication_Checker::decode_indices_from_index(int index, int & ai, int & aj, int & bi, int & bj)
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

inline int Matrix_Multiplication_Checker::get_vector_index(int a, int b)
{
    int result = a;
    result *= element_count;
    result += b;
    return result;
}

inline int Matrix_Multiplication_Checker::get_element_index(int i, int j)
{
    return ((i*length) + j);
}

inline int Matrix_Multiplication_Checker::get_m_index(int i, int j)
{
    return ((i*m_length) + j);
}

void Matrix_Multiplication_Checker::calculate_r_vectors()
{
    for (int i = 0; i < length; ++i)
        for (int j = 0; j < length; ++j)
        {
            int index = get_element_index(i, j);
            r_vectors[index].reset();
            for (int k = 0; k < length; ++k)
            {
                r_vectors[index][get_vector_index(i, k, k, j)] = 1;
            }
        }
}

void Matrix_Multiplication_Checker::calculate_m_vectors()
{
    for (int i = 1; i < power(2,element_count); ++i)
        for (int j = 1; j < power(2,element_count); ++j)
        {
            Bitvector_4 av(i);
            Bitvector_4 bv(j);
            int index = get_m_index(i-1, j-1);
            m_vectors[index].reset();
            //cout << i << " i= " << av << "\n" << j << " j= " << bv << "\n";
            for (int k = 0; k < element_count; ++k)
                for (int l = 0; l < element_count; ++l)
                {
                    if (av[k])
                        if (bv[l])
                            m_vectors[index].set(get_vector_index(k, l));
                }
        }
}

bool Matrix_Multiplication_Checker::check_for_good_vectors()
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
				vector<Bitvector_16> n_vectors;
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
				continue;
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

void Matrix_Multiplication_Checker::output_vector(Bitvector_16 v)
{
    for (int i = 0; i < v.size(); ++i)
    {
        cout << v[i];
    }
}

void Matrix_Multiplication_Checker::output_vector_text(Bitvector_16 v)
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

