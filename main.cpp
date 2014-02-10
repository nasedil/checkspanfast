#include <iostream>
#include <bitset>
#include <vector>
#include <fstream>
#include <ctime>
using namespace std;

const int MATRIX_SIZE = 2;

typedef bitset<MATRIX_SIZE*MATRIX_SIZE*MATRIX_SIZE*MATRIX_SIZE> Multiplication_vector;
typedef bitset<MATRIX_SIZE*MATRIX_SIZE> Matrix_vector;
typedef bitset<MATRIX_SIZE*MATRIX_SIZE*MATRIX_SIZE-1> Linear_vector;

// for positive b
int power(int a, int b);
// finds if there is at least one solution for system of linear equations in a vector form
bool gauss_solve(vector<Multiplication_vector> a, Multiplication_vector b);
// finds if there is at least one solution for system of linear equations in a vector form
bool binary_solve(vector<Multiplication_vector> a, Multiplication_vector b);
bool binary_solve_recursive(vector<Multiplication_vector> a, Multiplication_vector b);
bool binary_solve_explore(vector<Multiplication_vector> * a, Multiplication_vector * b, int depth, Multiplication_vector * current);

class Matrix_Multiplication_Checker
{
    public:
    int length;
    int element_count;
    Multiplication_vector * r_vectors;
    Multiplication_vector * m_vectors;
    int m_count;
    int m_length;
    int f_count;
    vector<Multiplication_vector> good_vectors;
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
    void output_vector(Multiplication_vector v);
    // outputs vector to screen in letters
    void output_vector_text(Multiplication_vector v);
    // finds if there are good vectors (what we want to find)
    bool check_for_good_vectors();
};

int main()
{
    // generate all possiple multiplications (> 200 vectors)

    Matrix_Multiplication_Checker * checker = new Matrix_Multiplication_Checker(MATRIX_SIZE);

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
    r_vectors = new Multiplication_vector[length*length];
    calculate_r_vectors();
    m_length = power(2,element_count)-1;
    m_count = m_length*m_length;
    m_vectors = new Multiplication_vector[m_count];
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
            Matrix_vector av(i);
            Matrix_vector bv(j);
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
  clock_t timestamp = clock();
    // for 2x2 matrix multiplication
    for (int c1 = 0; c1 < m_count-2; ++c1)
	{
		double elapsed_secs = double(clock() - timestamp) / CLOCKS_PER_SEC;
		timestamp = clock();
		cout << "time: " << elapsed_secs << "\n";
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
				vector<Multiplication_vector> n_vectors;
				// 9  37  62  109  128  146  165
				n_vectors.push_back(m_vectors[c1]);
				n_vectors.push_back(m_vectors[c2]);
				n_vectors.push_back(m_vectors[c3]);
				for (int i = 0; i < element_count; ++i)
					n_vectors.push_back(r_vectors[i]);
				int n_index = 0;
				for (int i = n_index; i < m_count; ++i)
				{
					if (binary_solve_recursive(n_vectors, m_vectors[i]))
					{
						if (good_vectors.size() == 0)
						{
							good_vectors.push_back(m_vectors[i]);
							good_vectors_indexes.push_back(i);
						}
						else
						{
							if (!binary_solve_recursive(good_vectors,m_vectors[i]))
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
				}
			}
		}
	}
}

void Matrix_Multiplication_Checker::output_vector(Multiplication_vector v)
{
    for (int i = 0; i < v.size(); ++i)
    {
        cout << v[i];
    }
}

void Matrix_Multiplication_Checker::output_vector_text(Multiplication_vector v)
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

int power(int a, int b)
{
    int result = a;
    while (--b)
        result *= a;
    return result;
}

bool gauss_solve(vector<Multiplication_vector> a, Multiplication_vector b)
{
    for (int i = 0; i < a.size(); ++i) // for all columns (7)
    {
        // first we exchange rows if needed, so that a[i][i] = 1
        if (!a[i][i])
        {
            int k = i+1;
            while (!a[i][k])
                ++k;
            bool tmp = b[k];
            b[k] = b[i];
            b[i] = tmp;
            for (int l = i; l < a.size(); ++l)
            {
                tmp = a[l][k];
                a[l][k] = a[l][i];
                a[l][i] = tmp;
            }
        }
        // second, we make all coefficients under a[i][i] equal to 0
        for (int j = i+1; j < b.size(); ++j) // for all rows (16)
        {
            if (a[i][j])
            {
                for (int k = i; k < a.size(); ++k)
                {
                    a[k][j] = (a[k][j] != a[k][i]);
                }
                b[j] = (b[j] != b[i]);
            }
        }
    }
    // we chech all i from 7 to 15 for zero coefficients
    for (int i = a.size(); i < b.size(); ++i)
    {
        if (b[i])
            return false;
    }
    return true;
}

bool binary_solve(vector<Multiplication_vector> a, Multiplication_vector b)
{
    uint_least64_t counter;
    uint_least64_t limit = power(2,MATRIX_SIZE*MATRIX_SIZE*MATRIX_SIZE-1);
    for (counter = 1; counter < limit; ++counter)
    {
        Linear_vector x(counter);
        Multiplication_vector r(0);
        bool is_good = true;
        // we multiply a by x
        for (int i = 0; i < b.size(); ++i)
        {
            for (int j = 0; j < a.size(); ++j)
            {
                if (x[j])
                    r[i] = (r[i] != a[j][i]);
            }
            if (r[i] != b[i])
            {
                is_good = false;
                break;
            }
        }
        if (is_good)
            return true;
    }
    return false;
}

bool binary_solve_recursive(vector<Multiplication_vector> a, Multiplication_vector b)
{
    uint_least64_t depth = a.size()-1;
    Multiplication_vector * current = new Multiplication_vector();
    current->reset();
	bool result = binary_solve_explore(&a, &b, depth, current);
	//cout << "\n!!!!";
    delete current;
    //cout << "\n!!!!";
    return result;
}

bool binary_solve_explore(vector<Multiplication_vector> * a, Multiplication_vector * b, int depth, Multiplication_vector * current)
{
	//cout << "\n" << depth;
	if (depth == -1)
	{
		//cout << "\nGo up";
		for (int k = 0; k < 16; ++k)
			if (current->test(k) != b->test(k))
				return false;
		return true;
		//return ((*current) == (*b));
	}
	else
	{
		//cout << "\nGo left";
		bool result = binary_solve_explore(a, b, depth-1, current);
		if (result)
			return true;
		//cout << "\nGo right";
		Multiplication_vector * current_add = new Multiplication_vector();
		for (int k = 0; k < 16; ++k)
			(*current_add)[k] = current->test(k) != (*a)[depth][k];
		result = binary_solve_explore(a, b, depth-1, current_add);
		delete current_add;
		//cout << "\nGo up";
		return result;
	}
}
