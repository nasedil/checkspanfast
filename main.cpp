#include <iostream>
#include <bitset>
#include <vector>
#include <fstream>
#include <ctime>
#include <algorithm>
#include "types.hpp"
#include "slae.hpp"
#include "Matrix_Multiplication_Checker.h"
using namespace std;

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


