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
#define VERY_DETAILED_OUTPUT /// output even more information
#define OUTPUT_SOLUTIONS_TO_FILE /// output solutions to a file
#define OUTPUT_STATISTICS

#include <iostream>

#include <omp.h>

#include "slae.hpp"
#include "cpchecker.hpp"

int main(int argc ,char** argv)
{
    if (argc > 1) {
        int p = atoi(argv[1]);	// processors
        if (omp_get_num_procs() < p) {
                std::cout << p << " processors are not available for this machine." << std::endl;
                return 0;
        } else {
            omp_set_num_threads(p);
        }
        #pragma omp parallel
        {
            //Cube_Product_Checker<2, 2, 16, 4>* checker = new Cube_Product_Checker<2, 2, 16, 4>;
            Cube_Product_Checker<2, 3, 512, 8>* checker = new Cube_Product_Checker<2, 3, 512, 8>;
            checker->init();
            if (checker->check_for_good_vectors_randomized()) {
                std::cout << "Found after " << checker->iteration_count << " iterations.";
            }
        }
    } else {
        std::cout << "Cube Product Checker (version pre-alpha)" << std::endl << std::endl
                  << "Usage:  cpc <threads>" << std::endl << std::endl
                  << "Starting default action..." << std::endl;

        Cube_Product_Checker<2, 2, 16, 4>* checker = new Cube_Product_Checker<2, 2, 16, 4>;
        //Cube_Product_Checker<2, 3, 512, 8>* checker = new Cube_Product_Checker<2, 3, 512, 8>;
        checker->init();
        checker->check_for_good_vectors();
        checker->save_results("./solutions-2x2.txt");
        bool all_right = true;
        Solution_Properties sp;
        Solution_Properties best_solution;
        best_solution.operation_count = 1000000;
        for (auto s: checker->solutions) {
            bool is_real = checker->check_solution(s, sp);
            if (best_solution.operation_count > sp.operation_count) {
                best_solution = sp;
            }
            if (!is_real) {
                all_right = false;
            }
        }
        if (all_right) {
            std::cout << "All solutions are checked and right." << std::endl;
        } else {
            std::cout << "All solutions are checked and there are problems!" << std::endl;
        }
        checker->save_solution_properties(best_solution, "./solution-2x2-best.txt");
        //checker->check_for_good_vectors_randomized();
        //checker->save_random_samples(1000000, "2test1.txt");
        //checker->read_samples_and_check("2test1.txt", "2test1-res.txt");
        delete checker;
    }

    return 0;
}
