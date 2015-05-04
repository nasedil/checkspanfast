/**
 * @file main.hpp:  Cube Product Checker main function.
 *
 * The program searches for Strassen-like algorithms on
 * 2- and 3-dimensional arrays with elements from GF(2).
 *
 * @author Eugene Petkevich
 * @version pre-alpha
 */

//#define VERBOSE_OUTPUT /// output verbose information
//#define VERY_DETAILED_OUTPUT /// output even more information
//#define OUTPUT_SOLUTIONS_TO_FILE /// output solutions to a file
//#define OUTPUT_STATISTICS
#define USE_CACHE

#include <iostream>
#include <string>

#include <omp.h>

#include "slae.hpp"
#include "cpchecker.hpp"
#include "utils.hpp"

using namespace std;

//=============================================================================

void temp(); /// tests during development process
void test(int n, int d, int t, int limit);  /// test different algorithms and write results
void parallel(int p, int n, int d, int limit); /// find solutions in parallel
void parallel_2x2();
void parallel_3x3(int limit);
void parallel_2x2x2(int limit);

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
        cout << "Cube Product Checker (version pre-alpha)" << endl << endl
                  << "Usage:  cpc <choice> <options>" << endl << endl
                  << "Choices are:" << endl
                  << "0  : test" << endl
                  << "1  : gather statistics" << endl
                  << "2  : run in parallel" << endl
                  << "3  : run distributed" << endl
                  << "4  : check solutions" << endl << endl;
    }
    return 0;
}

//=============================================================================
//=============================================================================

void temp()
{
    Cube_Product_Checker<2, 2, 16, 4>* checker = new Cube_Product_Checker<2, 2, 16, 4>;
    checker->init(7, "");
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
        checker->save_solution_properties(sp, string("./solution-2x2-" + to_string(i) + ".txt").c_str());
        if (!is_real) {
            all_right = false;
        }
        ++i;
    }
    if (all_right) {
        cout << "All solutions are checked and right." << endl;
    } else {
        cout << "All solutions are checked and there are problems!" << endl;
    }
    checker->save_solution_properties(best_solution, "./solution-2x2-best.txt");
    delete checker;
}

//=============================================================================

void parallel(int p, int n, int d, int limit)
{
    if (omp_get_num_procs() < p) {
            cout << p << " processors are not available for this machine." << endl;
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
            parallel_2x2x2(limit);
        }
    }
}

//=============================================================================

void parallel_2x2()
{
    #pragma omp parallel
    {
        Cube_Product_Checker<2, 2, 16, 4>* checker = new Cube_Product_Checker<2, 2, 16, 4>;
        checker->init(7, "");
        if (checker->check_for_good_vectors_randomized()) {
            cout << "Found after " << checker->iteration_count << " iterations.";
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
        //int thread_count = omp_get_num_threads();
        Cube_Product_Checker<3, 2, 81, 9>* checker;
        //if (omp_get_thread_num() == 0) {
        //    checker = master_checker;
        //} else {
            checker = new Cube_Product_Checker<3, 2, 81, 9>;
        //    checker->init(master_checker);
        checker->init(23, "");
        //}
        if (checker->solve_hill_climbing(limit)) {
            cout << "Found after " << checker->restarts << " restarts and " << checker->checked_sets_count << " checks.";
        }
        delete checker;
    }
}

void parallel_2x2x2(int limit)
{
    Cube_Product_Checker<2, 3, 512, 8>* master_checker = new Cube_Product_Checker<2, 3, 512, 8>;
    master_checker->init(15, "");
    //master_checker->write_m_vectors("mvectors.dat");

    #pragma omp parallel
    {
        int thread_number = omp_get_thread_num();
        Cube_Product_Checker<2, 3, 512, 8>* checker;
        if (thread_number == 0) {
            checker = master_checker;
        } else {
            checker = new Cube_Product_Checker<2, 3, 512, 8>;
            checker->init(*master_checker);
        }
        checker->thread_number = thread_number;
        if (checker->check_for_good_vectors_randomized()) {
            *(checker->stop_signal) = true;
            cout << "Found after " << checker->iteration_count << " iterations.";
        }
        // dump statistics
        if (thread_number > 0) {
            delete checker;
        }
    }
    #pragma omp barrier
    delete master_checker;
}

void test(int n, int d, int t, int limit)
{
    if ((n != 2) || (d != 2)) {
        cout << "The suite not ready yet.";
        return;
    }
    Cube_Product_Checker<2, 2, 16, 4>* checker = new Cube_Product_Checker<2, 2, 16, 4>;
    checker->init(7, "");
    Solution_Properties sp;
    Timewatch tw;
    int times = t;
    int checked_sets_count;
    int lin_dependent_sets;
    int iteration_count;
    int restarts;
    cout << "Start testing different search algorithms..." << endl;
    //----- test full space
    tw.watch();
    checker->check_for_good_vectors();
    cout << "Full search: " << tw.watch() << " s / "
              << checker->checked_sets_count << " checked sets" << endl;
    cout << "    " << checker->lin_dependent_sets << " linearly dependent sets" << endl;
    cout << "    " << checker->bit_check_hits << " bit check hits" << endl;
    cout << "    " << checker->gaussian_eliminations << " gaussian eliminations" << endl;
#ifdef USE_CACHE
    cout << "    " << checker->cache_hits << " cache hits" << endl;
#endif // USE_CACHE
    //----- test random search
    tw.watch();
    checked_sets_count = 0;
    lin_dependent_sets = 0;
    for (int i = 0; i < times; ++i) {
        checker->check_for_good_vectors_randomized();
        checked_sets_count += checker->checked_sets_count;
        lin_dependent_sets += checker->lin_dependent_sets;
#ifdef OUTPUT_STATISTICS
        if (!checker->check_solution(checker->good_vectors_indexes, sp)) {
            cout << "Found solution is not correct!!!" << endl;
            cout << "  Result is [ ";
            for (int i: checker->good_vectors_indexes) {
                cout << i << " ";
            }
            cout << "]" << endl;
        }
#endif // OUTPUT_STATISTICS
    }
    cout << "Random search: " << tw.watch()/times << " s / "
              << 1.0*checked_sets_count/times << " checked sets" << endl;
    cout << "    " << 1.0*lin_dependent_sets/times << " linearly dependent sets" << endl;
#ifdef USE_CACHE
    cout << "    " << checker->cache_hits << " cache hits" << endl;
#endif // USE_CACHE
    //----- test hill climbing
    tw.watch();
    checked_sets_count = 0;
    lin_dependent_sets = 0;
    iteration_count = 0;
    restarts = 0;
    for (int i = 0; i < times; ++i) {
        checker->solve_hill_climbing(limit);
        checked_sets_count += checker->checked_sets_count;
        lin_dependent_sets += checker->lin_dependent_sets;
        restarts += checker->restarts;
        iteration_count += checker->iteration_count;
#ifdef OUTPUT_STATISTICS
        if (!checker->check_solution(checker->good_vectors_indexes, sp)) {
            cout << "Found solution is not correct!!!" << endl;
            cout << "  Result is [ ";
            for (int i: checker->good_vectors_indexes) {
                cout << i << " ";
            }
            cout << "]" << endl;
        }
#endif // OUTPUT_STATISTICS
    }
    cout << "Hill climbing with restarts: " << tw.watch()/times << " s / "
              << 1.0*checked_sets_count/times << " checked sets / "
              << 1.0*restarts/times << " restarts / "
              << 1.0*iteration_count/times << " iterations" << endl;
    cout << "    " << 1.0*lin_dependent_sets/times << " linearly dependent sets" << endl;
#ifdef USE_CACHE
    cout << "    " << checker->cache_hits << " cache hits" << endl;
#endif // USE_CACHE
    //----- end
    cout << "Testing finished." << endl;
    delete checker;
}
