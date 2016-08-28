/*
 * NOTE: You must compile with the -lgmp, -lgmpxx, and -std=c++11 options or the linkage will fail.
 * Performs a series of calculations involving large numbers and wily functions.
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
#include <boost/multiprecision/gmp.hpp>
    using namespace boost::multiprecision;

#include "hpmath.hpp"

//////////////////////
///// prototypes /////
//////////////////////

void        console_io_exponential(void);
void        console_io_factorial(void);
void        console_io_partition(void);
void        console_io_tetration(void);
void        display_menu(void);
template <typename T>
bool        getInput(std::istream&, T&);
template <typename T>
bool        getRangedInput(std::istream&, T&, const T, const T);
MenuPick    read_choice(std::istream&);

//////////////////
///// main() /////
//////////////////

int main(void) {
    MenuPick selection;
    do {
        display_menu();
        selection = read_choice(cin);
        if (selection == tetration) {
            console_io_tetration();
        }
        else if (selection == factorials) {
            console_io_factorial();
        }
        else if (selection == exponentiate) {
            console_io_exponential();
        }
        else if (selection == partition) {
            console_io_partition();
        }
        else {
            selection = quit;
        }
    } while (selection != quit);

    return 0;
}

///////////////////////////////////
///// Free-floating functions /////
///////////////////////////////////

void console_io_exponential(void) {
    mpf_float_1000 x = 0.0;
    //unsigned int terms = 0;
    mpf_float_1000 result = 0.0;
    cout << "\nEnter the number to exponentiate: ";
    bool good;
    do {
        good = getInput(cin, x);
        if (!good) {
            cout << "Please enter a number: ";
        }
    } while (!good);
    result = exp(x);
    cout << "The result of exp(" << x << ") is " << std::setprecision(50)
         << result << "\n\n";
}

/*
 * get and display a tetration to the standard console i/o
 */
void console_io_factorial(void) {
    mpz_int n = 0;
    mpz_int result = 0;
    cout << "\nPlease enter an integer: ";
    bool good;
    do {
        good = getInput(cin, n);
        if (!good || n < 0) {
            cout << "Please enter a nonnegative integer: ";
        }
    } while (!good || n < 0);
    result = factorial(n);
    cout << "The result of " << n << "! is " << result << "\n\n";
}

/*
 * Performs the necessary I/O and calculations for the electrode
 */
void console_io_partition(void) {
    // data management
    unsigned int steps = 0; // the number of samples to take
    partition_fxn_sample* sample = nullptr; // sample array
    std::ofstream file; // the file to save to
    std::string filename; // the name of the file to save the results to
    progressBar<unsigned int> bar(100); // progress bar manipulator
    // variables
    unsigned short int states; // the number of states in the partition function
    mpf_float_1000* E = nullptr; // array of state energies
    mpf_float_1000 T = 0.0; // (K) Kelvin temperature (LCV)
    double T_max = 0.0; // (K) maximum Kelvin temperature
    double T_min = 0.0; // (K) minimum Kelvin temperature
    double step_size = 0.0; // (K) temperature step size
    // constants
    const mpf_float_1000 k_B = 8.6173324e-5; // (eV/K) Boltzmann constant

    // whether the stream input was acceptable
    bool good;
    cout << "\nEnter a filename to save the results (CSV format, will be overwritten): ";
    std::getline(cin, filename);
    cout << "How many states does the partition function have? ";
    getInput(cin, states);
    cout << "What is the maximum temperature to calculate (in K)?: ";
    do {
        good = getRangedInput(cin, T_max, 1e-100, std::numeric_limits<double>::max());
        if (!good) {
            cout << "Please enter a finite, positive temperature: ";
        }
    } while (!good);
    cout << "What is the minimum temperature to calculate (in K)?: ";
    do {
        good = getRangedInput(cin, T_min, 1e-100, T_max);
        if (!good || T_min == T_max) {
            cout << "Please enter a finite, positive temperature less than the maximum temperature: ";
        }
    } while (!good);
    cout << "How many Kelvins should the program step for each sample? ";
    do {
        good = getRangedInput(cin, step_size, 1e-100, T_max - T_min);
        if (!good) {
            cout << "Please enter a finite, positive value less than the temperature range: ";
        }
    } while (!good);

    // calculate the number of samples
    steps = static_cast<unsigned int>((T_max - T_min) / step_size);

    // allocate the energy array
    E = new mpf_float_1000[states];
    // allocate the sample array
    sample = new partition_fxn_sample[steps];

    // get the energies
    for (unsigned int i = 0; i < states; i++) {
        cout << "Enter the energy of the " << i+1
             << ((i+1 % 10 == 1 && i+1 % 100 != 11) ? "st" :
                 ((i+1 % 10 == 2 && i+1 % 100 != 12) ? "nd" :
                  ((i+1 % 10 == 3 && i+1 % 100 != 13) ? "rd" : "th")))
             << " state in eV: ";
        do {
            good = getInput(cin, E[i]);
            if (!good) {
                cout << "Please enter a number: ";
            }
        } while (!good);
    }

    cout << "Please wait . . .\n";
    // display the progress bar
    bar.initialize(cout, steps);
    
    // initialize the array
    for (unsigned int i = 0; i < steps; i++) {
        sample[i].P = new mpf_float_1000[states];
        sample[i].n = states;
        for (unsigned int j = 0; j < states; j++) {
            sample[i].P[j] = 0.0;
        }
    }

    // start at the minimum temp
    T = T_min;
    // acquire each sample
    for(unsigned int i = 0; i < steps; i++) {
        if (i != 0) {
            T += step_size;
        }
        sample[i].tau = k_B * T;
        // calculate the partition function; Z(tau) == \Sum_{j=0}^{\Infinity} \exp{-E / \tau}
        for (unsigned int j = 0; j < states; j++) {
            sample[i].P[j] = exp(-E[j] / sample[i].tau);
            sample[i].Z += sample[i].P[j];
        }
        // divide by Z to get the probability for each state
        for (unsigned int j = 0; j < states; j++) {
            sample[i].P[j] /= sample[i].Z;
        }
        
        // let the user know the program hasn't died yet
        bar.increment(i);
    }
    cout << "\nDone! Saving . . .\n\n";

    file.open(filename.c_str(), std::ofstream::out);
    if (file.is_open()) {
        file << std::setprecision(16) << "All energies are in eV\n\ntau,Z(tau)";
        // output the heading
        for (unsigned int i = 0; i < states; i++) {
            file << ",P_" << i+1;
        }
        file << '\n'; // start on the next row
        // output the data for each sample
        for (unsigned int i = 0; i < steps; i++) {
            file << sample[i].tau << ',' << sample[i].Z; // save tau and Z(tau)
            for (unsigned int j = 0; j < states; j++) {
                file << ',' << sample[i].P[j]; // output the probabilities
            }
            file << '\n'; // next row
        }

        file.close();
    }
    else {
        cout << "File error!\n";
    }
    
    // deallocate everything
    delete [] sample; // the destructors take care of the objects' pointers
    sample = nullptr;
}

/*
 * get and display a tetration to the standard console i/o
 */
void console_io_tetration(void) {
    mpf_float_1000 base = 0.0;
    int hyperpower = 0;
    mpf_float_1000 tetration_result = 0.0;
    bool good;
    cout << "\nPlease enter a base (a): ";
    do {
        good = getInput(cin, base);
        if (!good) {
            cout << "Please enter a number: ";
        }
    } while (!good);
    cout << "Please enter a hyperpower (n): ";
    do {
        good = getInput(cin, hyperpower);
        if (!good) {
            cout << "Please enter an integer: ";
        }
    } while (!good);
    tetration_result = tetrate(base, hyperpower);
    cout << "The result of " << base << "^^" << hyperpower
         << " is " << tetration_result << "\n\n";
}

/*
 * displays the selection menu
 */
void display_menu(void) {
    cout << "Please select an option by typing its number:\n"
         << "1. Tetrate a real number\n"
         << "2. Calculate a factorial\n"
         << "3. Exponentiate a real number\n"
         << "4. Evaluate the double-polymer electrode's partition function\n"
         << "5. Exit program\n";
}

/*
 * safely extract validated input from a stream
 * @param input         the stream to read from
 * @param value         the variable to save the data to
 * @return              whether or not the stream extraction was successful
 */
template <typename T>
bool getInput(std::istream& input, T& value) {
    bool success = true;
    input >> value;

    if (!input.good()) {
        success = false;
        // clear the flags
        input.clear();
        // empty the stream
        input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    else {
        // empty the stream
        input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    return success;
}

/*
 * safely extract input from a stream and ensure that it falls within a certain range
 * @param input         the stream to read from
 * @param value         the variable to save to
 * @param min           the minimum acceptable value (inclusive)
 * @param max           the maximum acceptable value (inclusive)
 */
template <typename T>
bool getRangedInput(std::istream& input, T& value, const T min, const T max) {
    bool success = getInput<T>(input, value);
    if (success) {
        if (value < min || value > max) {
            success = false;
        }
    }

    return success;
}

/*
 * get the selection from the user (or a file, in principle)
 * @param input         the input stream to read
 * @return              the selection
 */
MenuPick read_choice(std::istream& input) {
    MenuPick sel;
    int entered;
    bool good;
    do {
        good = getRangedInput(cin, entered, 1, 5);
        if (!good) {
            cout << "Please enter a number between 1 and 5: ";
        }
    } while (!good);
    if (entered == 1) {
        sel = tetration;
    }
    else if (entered == 2) {
        sel = factorials;
    }
    else if (entered == 3) {
        sel = exponentiate;
    }
    else if (entered == 4) {
        sel = partition;
    }
    else {
        sel = quit;
    }

    return sel;
}
