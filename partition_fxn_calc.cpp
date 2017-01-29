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

enum MenuChoice {
    varyTemp,
    varyVoltage,
    varyMagnet,
    quit
};

template <typename Num>
void sweepTemperature(Thermodynamics::SystemManager<Num>&);
template <typename Num>
void sweepElectricField(Thermodynamics::SystemManager<Num>&);

//////////////////
///// main() /////
//////////////////

int main(void) {
    Thermodynamics::SystemManager<mpfr_float_1000> system;
    system.params.acquire("config.cfg");

    sweepTemperature(system);

    cout << "\nSaving...\n";

    bool success;
    unsigned short tries = 3;
    do {
        success = system.save_to_disk(system.params.filename);
        if (!success) {
            cout << "The file could not be saved. Please enter a different file name: ";
            cin  >> system.params.filename;
        }
        tries--;
    } while (!success && (tries > 0));

    return 0;
}

///////////////////////////
///// other functions /////
///////////////////////////

/**
 *
 */
template <typename Num>
void sweepTemperature(Thermodynamics::SystemManager<Num>& system) {
    progressBar<unsigned int> pbar(80);
    Num T_min, T_max, T_step, T_current;

    cout << "What is the minimum temperature to calculate? ";
    getRangedInput(cin, T_min, static_cast<Num>(1e-100), static_cast<Num>(1e100));
    cout << "What is the maximum temperature to calculate? ";
    getRangedInput(cin, T_max, static_cast<Num>(1e-100), static_cast<Num>(1e100));
    cout << "What should the temperature step size be? ";
    getRangedInput(cin, T_step, static_cast<Num>(1e-100), static_cast<Num>(1e100));

    const unsigned short int n_samp = static_cast<unsigned short int>(static_cast<Num>((T_max - T_min) / T_step));
    system.initialize(n_samp);

    // just set total potential to zero for now
    for (unsigned short i = 0; i < system.params.states(); i++) {
        system.params.set_mu(i, 0.0);
    }
    
    std::cout << "Please wait . . .\n";
    pbar.initialize(cout, n_samp);

    T_current = T_min - T_step;
    for (unsigned short i = 0; i < n_samp; i++) {
        T_current += T_step;
        system.params.set_T(T_current);
        system.sample[i].calculate(system.params);
        pbar.increment(i);
    }

    pbar.end();
}

/**
 *
 */
template <typename Num>
void sweepElectricField(Thermodynamics::SystemManager<Num>& system) {
    
}
