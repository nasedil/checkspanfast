#include <iostream>
#include <bitset>
#include <vector>
#include <fstream>
#include <ctime>
#include <algorithm>
#include "slae.hpp"
#include "mmchecker.hpp"
using namespace std;

int main()
{
    Matrix_Multiplication_Checker<2, 2, 16, 4> * checker = new Matrix_Multiplication_Checker<2, 2, 16, 4>;
    //Matrix_Multiplication_Checker<2, 3, 512, 8> * checker = new Matrix_Multiplication_Checker<2, 3, 512, 8>;
    checker->check_for_good_vectors();
    //checker->check_for_good_vectors_randomized();

    delete checker;

    return 0;
}


