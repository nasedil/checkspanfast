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
    Matrix_Multiplication_Checker * checker = new Matrix_Multiplication_Checker(2);
    checker->check_for_good_vectors();

    delete checker;

    return 0;
}


