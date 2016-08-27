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

/////////////////////////////
///// Enums and Objects /////
/////////////////////////////

enum MenuPick {
    tetration,
    factorials,
    exponentiate,
    partition,
    quit
};

class partition_fxn_sample {
  public:
    ~partition_fxn_sample() {
        // clean up
        delete [] this->P;
        this->P = nullptr;
        this->n = 0;
    }
    mpf_float_1000 tau; // fundamental temperature
    mpf_float_1000 Z; // partition function at tau
    mpf_float_1000* P; // owning pointer to state probability array
    unsigned int n; // size of state probability array
    partition_fxn_sample(const unsigned int numstates) {
        // allocate the probability state array
        this->P = new mpf_float_1000[numstates];
        this->n = numstates;
    }
    partition_fxn_sample(void) {
        // empty
        this->tau = this->Z = 0.0;
        this->P = nullptr;
        this->n = 0;
    }
};

template <typename Num>
class progressBar {
    Num full;
    Num current;
    std::ostream* stream;
    unsigned int width;
  public:
    progressBar(void) {
        this->full = this->current = 0;
        this->stream = nullptr;
        this->width = 80;
    }
    progressBar(unsigned int w) {
        this->full = this->current = 0;
        this->stream = nullptr;
        this->width = w;
    }
    void initialize(std::ostream& output, const Num total) {
        // set the 100% mark
        this->full = total;
        // save the pointer to the stream
        this->stream = &output;
        // print out the initial, empty indicator
        *(this->stream) << '[';
        for (unsigned int i = 0; i < (this->width-2); i++) {
            *(this->stream) << ' ';
        }
        // cap off the progress bar and return to the beginning of the line
        *(this->stream) << "]\r[|";
        // flush the buffer to make sure everything prints
        this->stream->flush();
    }
    void increment(const Num now) {
        Num frac = static_cast<Num>((this->width) * static_cast<double>(now) / this->full);
        Num diff = frac - this->current;
        for (Num i = 0; i < diff; i++) {
            // increment with a pipe
            *(this->stream) << '|';
            // ensure it's printed
            this->stream->flush();
        }
        this->current = frac;
    }
    ~progressBar(void) {
        this->full = this->current = 0;
        this->stream = nullptr;
    }
};

//////////////////////
///// prototypes /////
//////////////////////

void            console_io_exponential(void);
void            console_io_factorial(void);
void            console_io_partition(void);
void            console_io_tetration(void);
void            display_menu(void);
template <typename Numerical>
Numerical       exp(const Numerical);
template <typename Numerical>
Numerical       factorial(const Numerical);
void            printequals(std::ostream&, const unsigned int);
MenuPick        read_choice(std::istream&);
template <typename Numerical>
Numerical       tetrate(const Numerical, const unsigned int);

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

/////////////////////////////
///// Support Functions /////
/////////////////////////////

void console_io_exponential(void) {
    mpf_float_1000 x = 0.0;
    //unsigned int terms = 0;
    mpf_float_1000 result = 0.0;
    cout << "Enter the number to exponentiate: ";
    cin  >> x;
    result = exp(x);
    cout << "The result of exp(" << x << ") is " << std::setprecision(50) << result << "\n";
}

/*
 * get and display a tetration to the standard console i/o
 */
void console_io_factorial(void) {
    mpz_int n = 0;
    mpz_int result = 0;
    cout << "Please enter an integer: ";
    cin  >> n;
    result = factorial(n);
    cout << "The result of " << n << "! is " << result << "\n";
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
    unsigned int states; // the number of states in the partition function
    mpf_float_1000* E = nullptr; // array of state energies
    mpf_float_1000 T = 0.0; // (K) Kelvin temperature (LCV)
    double T_max = 0.0; // (K) maximum Kelvin temperature
    double T_min = 0.0; // (K) minimum Kelvin temperature
    double step_size = 0.0; // (K) temperature step size
    // constants
    const mpf_float_1000 k_B = 8.6173324e-5; // (eV/K) Boltzmann constant

    cout << "\nEnter a filename to save the results (CSV format, will be overwritten): ";
    std::getline(cin, filename);
    cout << "How many states does the partition function have? ";
    cin  >> states;
    cout << "What is the maximum temperature to calculate (in K)?: ";
    cin  >> T_max;
    cout << "What is the minimum temperature to calculate (in K)?: ";
    cin  >> T_min;
    cout << "How many Kelvins should the program step for each sample? ";
    cin  >> step_size;

    // don't allow absolute zero (mathematical singularity exists)
    if (T_min == 0.0) {
        T_min = 1e-10;
    }
    // don't allow negative temps
    if (T_min < 0.0) {
        T_min = 1e-10;
    }
    if (T_max < 0.0) {
        T_max = T_min + 10 * step_size;
    }
    // if the max and min are backwards, silently fix it
    if (T_min > T_max) {
        double temp = T_min;
        T_min = T_max;
        T_max = temp;
    }

    // calculate the number of samples
    steps = static_cast<unsigned int>((T_max - T_min) / step_size);

    // allocate the energy array
    E = new mpf_float_1000[states];
    // allocate the sample array
    sample = new partition_fxn_sample[steps];

    // get the energies
    for (unsigned int i = 0; i < states; i++) {
        cout << "Enter the energy of the " << i+1
             << (i+1 == 1 ? "st" : (i+1 % 10 == 2 ? "nd" : (i+1 % 10 == 3 ? "rd" : "th")))
             << " state in eV: ";
        cin  >> E[i];
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
    unsigned int hyperpower = 0;
    mpf_float_1000 tetration_result = 0.0;
    cout << "Please enter a base (a): ";
    cin  >> base;
    cout << "Please enter a hyperpower (n): ";
    cin  >> hyperpower;
    tetration_result = tetrate(base, hyperpower);
    cout << "The result of " << base << "^^" << hyperpower
         << " is " << tetration_result << "\n";
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
 * Raises Euler's number to a power (doesn't stop until all digits are filled!)
 * @param x             the number to raise e to
 * @return              the result
 */
template <typename Numerical>
Numerical exp(const Numerical x) {
    Numerical result = 1.0; // 0th term == 1
    Numerical prev = 0.0;
    unsigned long long int i = 1;
    // e^x == \Sum_{n==0}^{Inf} \frac{x^n}{n!}
    while (result != prev && i <= std::numeric_limits<unsigned long long int>::max()) {
        prev = result;
        Numerical numerator = 1.0;
        // x^i on index i without using double pow(mpf_float_1000, mpz_int)
        for (unsigned long long int j = 1; j <= i; j++) {
            numerator *= x;
        }
        result += numerator / factorial(i);
        i++;
    }

    return result;
}

/*
 * calculate a factorial of a nonnegative integer
 */
template <typename Numerical>
Numerical factorial(const Numerical n) {
    Numerical result = 1;
    // will just return 1 if the argument is < 2
    for (Numerical i = 2; i <= n; i++) {
        result *= i;
    }

    return result;
}

/*
 * print a number of equals signs to output
 * @param output        the output stream to write equals signs to
 * @param num           the numer of '=' to print
 */
void printequals(std::ostream& output, const unsigned int num) {
    for (unsigned int i = 0; i < num; i++) {
        output << '=';
    }
}

/*
 * get the selection from the user (or a file, in principle)
 * @param input         the input stream to read
 * @return              the selection
 */
MenuPick read_choice(std::istream& input) {
    MenuPick sel;
    int entered;
    input >> entered;
    input.ignore(1);
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

/*
 * A template for a tetration function which calculates a power tower for
 * positive hyperpowers and repeated roots for negative hyperpowers
 * @ param base         the number to be tetrated
 * @ param hyperpower   the height of the power/root tower
 * @ return             the result of the tetration
 */
template <typename Numerical>
Numerical tetrate(const Numerical base, const unsigned int hyperpower) {
    Numerical result = base;
    // repeated exponentiation for positive hyperpowers
    if (hyperpower > 0) {
        for (int i = 1; i < hyperpower; i++) {
            result = pow(base, result);
        }
    }
    // repeated roots for negative powers
    else if (hyperpower < 0 && base >= 0) {
        // decrement since the hyperpower is negative
        for (int i = -1; i > hyperpower; i--) {
            result = pow(base, 1/result);
        }
    }
    // return 1.00 for zero powers
    else {
        result = 1;
    }
    
    return result;
}
