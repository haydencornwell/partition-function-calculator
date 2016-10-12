/*
 * NOTE: You must compile with the -lmpfr and -std=c++11 options or the linkage will fail.
 * @author      Hayden Cornwell
 */

#include <iostream>
    using std::cin;
    using std::cout;
    using std::endl;
#include <iomanip>
#include <limits>
#include <string>
#include <fstream>
#include <cmath>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/cpp_int.hpp>
    using namespace boost::multiprecision;

#include "classes.hpp"

//////////////////
///// main() /////
//////////////////

int main(void) {
    mpfr_float_1000 T_current;
    progressBar<unsigned int> pbar(80);
    SystemManager<mpfr_float_1000> system;

    system.initialize("config.cfg");
    T_current = system.params.T_min() - system.params.T_step();
    pbar.initialize(cout, system.params.n_samp());

    for (unsigned int i = 0; i < system.params.n_samp(); i++) {
        T_current += system.params.T_step();
        system.sample[i].calculate(T_current, system.params);
        pbar.increment(i);
    }
    pbar.end();
    cout << "\nSaving...\n";

    bool success;
    do {
        success = system.save_to_disk(system.params.filename);
        if (!success) {
            cout << "The file could not be saved. Please enter a different file name: ";
            cin  >> system.params.filename;
        }
    } while (!success);

    return 0;
}
