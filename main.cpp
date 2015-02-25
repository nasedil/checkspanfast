//#define output2
//#define output1
//#define statistics

#include <iostream>

#include "slae.hpp"
#include "mmchecker.hpp"

int main(int argc , char * argv)
{
    if (argc > 1)
    {
        std::cout << "no arguments supported yet!";
    }
    else
    {
        std::cout << "Cube Product Checker (version pre-alpha)" << std::endl << std::endl
                  << "Usage:  ..." << std::endl << std::endl
                  << "Starting default action...";

        Matrix_Multiplication_Checker<2, 2, 16, 4> * checker = new Matrix_Multiplication_Checker<2, 2, 16, 4>;
        //Matrix_Multiplication_Checker<2, 3, 512, 8> * checker = new Matrix_Multiplication_Checker<2, 3, 512, 8>;
        checker->check_for_good_vectors();
        //checker->check_for_good_vectors_randomized();
        //checker->save_random_samples(1000000, "2test1.txt");
        //checker->read_samples_and_check("2test1.txt", "2test1-res.txt");
        delete checker;
    }

    return 0;
}


