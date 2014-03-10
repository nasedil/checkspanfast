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
    checker->check_for_good_vectors();

    delete checker;

    return 0;
}


