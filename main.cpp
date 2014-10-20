//#define output2
//#define output1

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

    checker->start_statistics("stats2d.txt");

    //checker->check_for_good_vectors();
    //checker->check_for_good_vectors_randomized();

    //checker->save_random_samples(1000000, "2test1.txt");
    checker->read_samples_and_check("2test1.txt", "2test1-res.txt");

	checker->end_statistics();

    delete checker;

    return 0;
}


