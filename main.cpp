#include <iostream>
#include <bitset>
using namespace std;

const int MATRIX_SIZE = 2;

typedef bitset<MATRIX_SIZE*MATRIX_SIZE*MATRIX_SIZE*MATRIX_SIZE> Multiplication_vector;

class Matrix_Multiplication_Checker
{
    public:
    int length;
    //int bit_for_length;
    int element_count;
    Multiplication_vector * r_vectors;

    Matrix_Multiplication_Checker(int length);
    ~Matrix_Multiplication_Checker();

    // returns linear index in a bit vector by 4 indexes in a matrices
    int get_vector_index(int ai, int aj, int bi, int bj);
    // returns linear index of an elementt in a matrix
    int get_element_index(int i, int j);
    // write result vectors to array
    void calculate_r_vectors();
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
    // calculate bit_for_length
    /*
    this->bit_for_length = 0;
    while (length/=2)
        ++(this->bit_for_length);
    */
    // calculate r_vectors
    r_vectors = new Multiplication_vector[length*length];
    calculate_r_vectors();
}

Matrix_Multiplication_Checker::~Matrix_Multiplication_Checker()
{
    delete [] r_vectors;
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

int Matrix_Multiplication_Checker::get_element_index(int i, int j)
{
    return ((i*length) + j);
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

void Matrix_Multiplication_Checker::output_vector(Multiplication_vector v)
{
    for (int i = 0; i < v.size(); ++i)
    {
        cout << v[i];
    }
}

