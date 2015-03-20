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
#define OUTPUT_SOLUTIONS_TO_FILE /// output solutions to a file
//#define OUTPUT_STATISTICS

#include <iostream>
#include <string>

#include <omp.h>

#include "slae.hpp"
#include "cpchecker.hpp"
#include "utils.hpp"

//=============================================================================

void temp(); /// tests during development process
void test(int n, int d, int t, int limit);  /// test different algorithms and write results
void parallel(int p, int n, int d, int limit); /// find solutions in parallel
void parallel_2x2();
void parallel_3x3(int limit);
void parallel_2x2x2();

//=============================================================================

int main(int argc ,char** argv)
{
    if (argc > 1) {
        int choice = atoi(argv[1]);
        switch (choice) {
            case 1: {
                int n = atoi(argv[2]);
                int d = atoi(argv[3]);
                int t = atoi(argv[4]);
                int limit = atoi(argv[5]);
                test(n, d, t, limit);
                break;
            }
            case 2: {
                int p = atoi(argv[2]);
                int n = atoi(argv[3]);
                int d = atoi(argv[4]);
                int limit = atoi(argv[5]);
                parallel(p, n, d, limit);
                break;
            }
            default: {
                temp();
                break;
            }
        }
    } else {
        std::cout << "Cube Product Checker (version pre-alpha)" << std::endl << std::endl
                  << "Usage:  cpc <choice> <options>" << std::endl << std::endl
                  << "Choices are:" << std::endl
                  << "0  : test" << std::endl
                  << "1  : gather statistics" << std::endl
                  << "2  : run in parallel" << std::endl
                  << "3  : run distributed" << std::endl
                  << "4  : check solutions" << std::endl << std::endl;
    }
    return 0;
}

//=============================================================================
//=============================================================================

void temp()
{
    Cube_Product_Checker<2, 2, 16, 4>* checker = new Cube_Product_Checker<2, 2, 16, 4>;
    checker->init(7);
    checker->check_for_good_vectors();
    checker->save_results("./solutions-2x2.txt");
    bool all_right = true;
    Solution_Properties sp;
    Solution_Properties best_solution;
    best_solution.operation_count = 1000000;
    int i = 1;
    for (auto s: checker->solutions) {
        bool is_real = checker->check_solution(s, sp);
        if (best_solution.operation_count > sp.operation_count) {
            best_solution = sp;
        }
        checker->save_solution_properties(sp, std::string("./solution-2x2-" + std::to_string(i) + ".txt").c_str());
        if (!is_real) {
            all_right = false;
        }
        ++i;
    }
    if (all_right) {
        std::cout << "All solutions are checked and right." << std::endl;
    } else {
        std::cout << "All solutions are checked and there are problems!" << std::endl;
    }
    checker->save_solution_properties(best_solution, "./solution-2x2-best.txt");
    delete checker;
}

//=============================================================================

void parallel(int p, int n, int d, int limit)
{
    if (omp_get_num_procs() < p) {
            std::cout << p << " processors are not available for this machine." << std::endl;
            return;
    } else {
        omp_set_num_threads(p);
    }
    if (d == 2) {
        if (n == 2) {
            parallel_2x2();
        } else if (n == 3) {
            parallel_3x3(limit);
        }
    } else if (d == 3) {
        if (n == 2) {
            parallel_2x2x2();
        }
    }
}

//=============================================================================

void parallel_2x2()
{
    #pragma omp parallel
    {
        Cube_Product_Checker<2, 2, 16, 4>* checker = new Cube_Product_Checker<2, 2, 16, 4>;
        checker->init(7);
        if (checker->check_for_good_vectors_randomized()) {
            std::cout << "Found after " << checker->iteration_count << " iterations.";
        }
        delete checker;
    };
}

//=============================================================================

void parallel_3x3(int limit)
{
    //Cube_Product_Checker<3, 2, 81, 9>* master_checker = new Cube_Product_Checker<3, 2, 81, 9>;
    //master_checker->init(26);
    #pragma omp parallel
    {
        int thread_count = omp_get_num_threads();
        Cube_Product_Checker<3, 2, 81, 9>* checker;
        //if (omp_get_thread_num() == 0) {
        //    checker = master_checker;
        //} else {
            checker = new Cube_Product_Checker<3, 2, 81, 9>;
        //    checker->init(master_checker);
        checker->init(23);
        //}
        checker->tw.thread_count = thread_count;
        if (checker->solve_hill_climbing(limit)) {
            std::cout << "Found after " << checker->restarts << " restarts and " << checker->checked_sets_count << " checks.";
        }
        delete checker;
    }
}

void parallel_2x2x2()
{
    Cube_Product_Checker<2, 3, 512, 8>* master_checker = new Cube_Product_Checker<2, 3, 512, 8>;
    master_checker->init(15);
    #pragma omp parallel
    {
        int thread_count = omp_get_num_threads();
        Cube_Product_Checker<2, 3, 512, 8>* checker;
        if (omp_get_thread_num() == 0) {
            checker = master_checker;
        } else {
            checker = new Cube_Product_Checker<2, 3, 512, 8>;
            checker->init(*master_checker);
        }
        checker->tw.thread_count = thread_count;
        if (checker->check_for_good_vectors_randomized()) {
            std::cout << "Found after " << checker->iteration_count << " iterations.";
        }
        delete checker;
    }
}

void test(int n, int d, int t, int limit)
{
    if ((n != 2) || (d != 2)) {
        std::cout << "The suite not ready yet.";
        return;
    }
    Cube_Product_Checker<2, 2, 16, 4>* checker = new Cube_Product_Checker<2, 2, 16, 4>;
    checker->init(7);
    Solution_Properties sp;
    Timewatch tw;
    int times = t;
    int checked_sets_count;
    int iteration_count;
    int restarts;
    std::cout << "Start testing different search algorithms..." << std::endl;
    //----- test full space
    //tw.watch();
    //checker->check_for_good_vectors();
    //std::cout << "Full search: " << tw.watch() << " s / "
    //          << checker->checked_sets_count << " checked sets" << std::endl;
    //----- test random search
    tw.watch();
    checked_sets_count = 0;
    for (int i = 0; i < times; ++i) {
        checker->check_for_good_vectors_randomized();
        checked_sets_count += checker->checked_sets_count;
#ifdef OUTPUT_STATISTICS
        if (!checker->check_solution(checker->good_vectors_indexes, sp)) {
            std::cout << "Found solution is not correct!!!" << std::endl;
            std::cout << "  Result is [ ";
            for (int i: checker->good_vectors_indexes) {
                std::cout << i << " ";
            }
            std::cout << "]" << std::endl;
        }
#endif // OUTPUT_STATISTICS
    }
    std::cout << "Random search: " << tw.watch()/times << " s / "
              << 1.0*checked_sets_count/times << " checked sets" << std::endl;
    //----- test hill climbing
    tw.watch();
    checked_sets_count = 0;
    iteration_count = 0;
    restarts = 0;
    for (int i = 0; i < times; ++i) {
        checker->solve_hill_climbing(limit);
        checked_sets_count += checker->checked_sets_count;
        restarts += checker->restarts;
        iteration_count += checker->iteration_count;
#ifdef OUTPUT_STATISTICS
        if (!checker->check_solution(checker->good_vectors_indexes, sp)) {
            std::cout << "Found solution is not correct!!!" << std::endl;
            std::cout << "  Result is [ ";
            for (int i: checker->good_vectors_indexes) {
                std::cout << i << " ";
            }
            std::cout << "]" << std::endl;
        }
#endif // OUTPUT_STATISTICS
    }
    std::cout << "Hill climbing with restarts: " << tw.watch()/times << " s / "
              << 1.0*checked_sets_count/times << " checked sets / "
              << 1.0*restarts/times << " restarts / "
              << 1.0*iteration_count/times << " iterations" << std::endl;
    //----- end
    std::cout << "Testing finished." << std::endl;
    delete checker;
}
