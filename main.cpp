/**
 * @file main.hpp:  Cube Product Checker main function.
 *
 * The program searches for Strassen-like algorithms on
 * 2- and 3-dimensional arrays with elements from GF(2).
 *
 * @author Eugene Petkevich
 * @version pre-alpha
 */

#define VERBOSE_OUTPUT /// output verbose information
//#define VERY_DETAILED_OUTPUT /// output even more information
//#define OUTPUT_SOLUTIONS_TO_FILE /// output solutions to a file
//#define OUTPUT_STATISTICS

#include <iostream>

#include "slae.hpp"
#include "mmchecker.hpp"

int main(int argc ,char ** argv)
{
    if (argc > 1)
    {
        std::cout << "no arguments supported yet!";
    }
    else
    {
        std::cout << "Cube Product Checker (version pre-alpha)" << std::endl << std::endl
                  << "Usage:  ..." << std::endl << std::endl
                  << "Starting default action..." << std::endl;

        Cube_Product_Checker<2, 2, 16, 4> * checker = new Cube_Product_Checker<2, 2, 16, 4>;
        //Cube_Product_Checker<2, 3, 512, 8> * checker = new Cube_Product_Checker<2, 3, 512, 8>;
        checker->check_for_good_vectors();
        //checker->check_for_good_vectors_randomized();
        //checker->save_random_samples(1000000, "2test1.txt");
        //checker->read_samples_and_check("2test1.txt", "2test1-res.txt");
        delete checker;
    }

    return 0;
}


