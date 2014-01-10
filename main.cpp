#include <iostream>
#include <bitset>
using namespace std;

const int MATRIX_SIZE = 2;

typedef bitset<MATRIX_SIZE*MATRIX_SIZE*MATRIX_SIZE*MATRIX_SIZE> Multiplication_vector;
typedef bitset<MATRIX_SIZE*MATRIX_SIZE> Matrix_vector;

// for positive b
int power(int a, int b);

class Matrix_Multiplication_Checker
{
    public:
    int length;
    int element_count;
    Multiplication_vector * r_vectors;
    Multiplication_vector * m_vectors;
    int m_count;
    int m_length;

    Matrix_Multiplication_Checker(int length);
    ~Matrix_Multiplication_Checker();

    // returns linear index in a bit vector by 4 indexes in matrices
    int get_vector_index(int ai, int aj, int bi, int bj);
    // returns linear index in a bit vector by 2 indexes in matrices
    int get_vector_index(int a, int b);
    // returns linear index of an elementt in a matrix
    int get_element_index(int i, int j);
    // returns index of vertar in m
    int get_m_index(int i, int j);
    // write result vectors to array
    void calculate_r_vectors();
    // write different multiplication vectors to array
    void calculate_m_vectors();
    // outputs vector to memory
    void output_vector(Multiplication_vector v);
};

int main()
{
    // generate all possiple multiplications (> 200 vectors)

    Matrix_Multiplication_Checker * checker = new Matrix_Multiplication_Checker(MATRIX_SIZE);

    for (int i = 0; i < checker->length; ++i)
        for (int j = 0; j < checker->length; ++j)
        {
            checker->output_vector(checker->r_vectors[checker->get_element_index(i, j)]);
            cout << "\n";
        }

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
}

Matrix_Multiplication_Checker::~Matrix_Multiplication_Checker()
{
    delete [] r_vectors;
    delete [] m_vectors;
}

int Matrix_Multiplication_Checker::get_vector_index(int ai, int aj, int bi, int bj)
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

int Matrix_Multiplication_Checker::get_vector_index(int a, int b)
{
    int result = a;
    result *= element_count;
    result += b;
    return result;
}

int Matrix_Multiplication_Checker::get_element_index(int i, int j)
{
    return ((i*length) + j);
}

int Matrix_Multiplication_Checker::get_m_index(int i, int j)
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
            cout << index << " : ";
            output_vector(m_vectors[index]);
            cout << "\n";
        }
}

void Matrix_Multiplication_Checker::output_vector(Multiplication_vector v)
{
    for (int i = 0; i < v.size(); ++i)
    {
        cout << v[i];
    }
}


int power(int a, int b)
{
    int result = a;
    while (--b)
        result *= a;
    return result;
}

