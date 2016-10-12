#ifndef CLASSES_HPP
    #define CLASSES_HPP

#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include <fstream>
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/cpp_int.hpp>
    using namespace boost::multiprecision;

#include "hpmath.hpp"
#include "templates.hpp"

/////////////////////
///// Constants /////
/////////////////////
namespace Constants {
    const mpfr_float_1000 k_B = 8.6173324e-5_mpr1k;             // (eV/K) Boltzmann constant
    const mpfr_float_1000 N_A = 6.022140857e-5_mpr1k;           // (per mol) Avogadro constant
    const mpfr_float_1000 V_abs = 4.44_mpr1k;                   // (V) absolute potential
    const mpfr_float_1000 mpfr_pi = boost::math::constants::pi<mpfr_float_1000>();
    const mpfr_float_1000 permittivity = 8.854187817e-12_mpr1k; // (F / m) electric permittivity
    const mpfr_float_1000 k_e = 1_mpr1k / (4_mpr1k * mpfr_pi * permittivity); // (N * m^2 / C^2) Coulomb's constant
}

///////////////////
///// Objects /////
///////////////////

template <typename Num>
class SystemParameters {
    unsigned short int n;   // size of the energy array
    Num* E;                 // (eV) owning pointer to energy array
    Num T_MIN;              // (K) minimum temperature
    Num T_MAX;              // (K) maximum temperature
    Num step_size;          // (K) temperature step size
  public:
    std::string filename;   // the name of the file to save to
    ~SystemParameters(void);
    SystemParameters(void);
    unsigned short int states(void) {return this->n;}
    unsigned int n_samp(void);
    Num T_min(void) {return this->T_MIN;}
    Num T_step(void) {return this->step_size;}
    Num energy(unsigned short int i) {return (i >= 0 && i < this->n ? this->E[i] : static_cast<Num>(0));}
    SystemParameters& acquire(const std::string);
};

template <typename Num>
class PartitionFunctionSample {
    unsigned short int states;  // size of state probability array
    Num TAU;                    // (eV) fundamental temperature
    Num PARTITION;              // partition function at tau
    Num* P;                     // owning pointer to state probability array
    Num TEMPERATURE;            // (K) temperature
  public:
    ~PartitionFunctionSample(void);
    PartitionFunctionSample(const unsigned int);
    PartitionFunctionSample(void);
    // accessors
    Num tau(void) {return this->TAU;}
    Num Z(void) {return this->PARTITION;}
    Num P_i(unsigned short int i) {return (i < this->states ? this->P[i] : 0);}
    Num T(void) {return this->TEMPERATURE;}
    // calculation
    void calculate(Num, SystemParameters<Num>&);
    // initialization in case the default constructor was used
    void initialize(unsigned short int);
};

template <typename Num>
class SystemManager {
  public:
    SystemParameters<Num> params;         // thermodynamic system parameters
    PartitionFunctionSample<Num>* sample; // owning pointer to sample array

    ~SystemManager(void);
    SystemManager(void);
    bool save_to_disk(std::string);
    void initialize(const std::string);
};

template <typename Num>
class progressBar {
    Num full;
    Num current;
    std::ostream* stream;
    unsigned int width;
  public:
    ~progressBar(void);
    progressBar(void);
    progressBar(unsigned int);
    void initialize(std::ostream&, const Num);
    void increment(const Num);
    void end(void);
};

///////////////////////////////////////
///// Member Function Definitions /////
///////////////////////////////////////

/* class SystemParameters */

template <typename Num>
SystemParameters<Num>::~SystemParameters(void) {
    delete [] this->E;
    this->n = 0;
}

template <typename Num>
SystemParameters<Num>::SystemParameters(void) {
    this->n = 0;
}

template <typename Num>
unsigned int SystemParameters<Num>::n_samp(void) {
    return static_cast<unsigned int>(static_cast<Num>((this->T_MAX - this->T_MIN) / this->step_size));
}

template <typename Num>
SystemParameters<Num>& SystemParameters<Num>::acquire(const std::string cfg_name) {
    std::ifstream config(cfg_name);
    std::string use_cfg_response;
    // try to read from a config file first
    if (config.is_open()) {
        cout << "\nConfiguration file found; use data? (y/n) ";
        std::getline(cin, use_cfg_response);
        if (static_cast<char>(tolower(use_cfg_response[0])) == 'y') {
            std::getline(config, this->filename, '\n');
            // handle Windows/DOS line endings when using std::getline
            if (this->filename.back() == '\r') {
                this->filename.pop_back(); // delete the last character, which is a carriage return
            }
            config >> this->T_MIN >> this->T_MAX >> this->n;
            if (this->n != 0) {
                E = new Num[this->n];
                for (unsigned short int i = 0; i < this->n; i++) {
                    config >> this->E[i];
                }
            }
            config.close();
        }
        else {
            config.close();
        }

        if (!config.good() || this->n == 0) {
            cout << "Error reading configuration file.\n\n";
            config.close();
        }
    }

    if (static_cast<char>(tolower(use_cfg_response[0])) == 'n' || !config.good()) {
        // file name
        cout << "\nEnter a filename to save the results (CSV format, will be overwritten): ";
        std::getline(cin, this->filename);
    
        // temperature
        cout << "What is the minimum temperature to calculate (in K)?: ";
        rangedGetterLoop(cin, cout, this->T_MIN, 1e-100_mpr1k, 1.416833e32_mpr1k,
                         "Please enter a finite, positive temperature in Kelvins: ");

        cout << "What is the maximum temperature to calculate (in K)?: ";
        rangedGetterLoop(cin, cout, this->T_MAX, this->T_MIN, 1.416833e32_mpr1k,
                         "Please enter a finite, positive temperature less than the minimum temperature: ");

        cout << "How many Kelvins should the program step for each sample? ";
        rangedGetterLoop(cin, cout, this->step_size, static_cast<Num>(1e-100),
                         static_cast<Num>(this->T_MAX - this->T_MIN),
                         "Please enter a finite, positive value less than the temperature range: ");

        cout << "How many states does the partition function have? ";
        rangedGetterLoop(cin, cout, this->n, static_cast<unsigned short>(0),
                         std::numeric_limits<unsigned short int>::max(),
                         "Please enter a positive integer: ");

        E = new Num[this->n];

        // get the energies
        for (unsigned short int i = 0; i < this->n; i++) {
            cout << "Enter the energy of the " << i+1
                 << ((i+1 % 10 == 1 && i+1 % 100 != 11) ? "st" :
                    ((i+1 % 10 == 2 && i+1 % 100 != 12) ? "nd" :
                    ((i+1 % 10 == 3 && i+1 % 100 != 13) ? "rd" : "th")))
                 << " state in eV: ";
            getterLoop(cin, cout, this->E[i], "Please enter a numerical value: ");
        }

        cout << "Please wait . . .\n";
    }

    // returns itself so this can be called inside of another function that uses this type.
    return *this;
}

/* class PartitionFunctionSample */

template <typename Num>
PartitionFunctionSample<Num>::~PartitionFunctionSample(void) {
    // clean up
    delete [] this->P;
    this->P = nullptr;
    this->states = 0;
}

template <typename Num>
PartitionFunctionSample<Num>::PartitionFunctionSample(void) {
    // empty
    this->TAU = this->PARTITION = 0.0;
    this->P = nullptr;
    this->states = 0;
}

template <typename Num>
PartitionFunctionSample<Num>::PartitionFunctionSample(const unsigned int numstates) {
    this->TAU = this->PARTITION = 0.0;
    // allocate the probability state array
    this->P = new mpfr_float_1000[numstates];
    this->states = numstates;
    for (unsigned int i = 0; i < this->states; i++) {
        this->P[i] = 0.0;
    }
}

/*
 * calculate the values at the given temperature and state energies
 * @param temp          the current temperature in Kelvins
 * @param E             a pointer to an array of state energy values in eV
 */
template <typename Num>
void PartitionFunctionSample<Num>::calculate(Num temp, SystemParameters<Num>& E) {
    this->TEMPERATURE = temp;
    this->TAU = Constants::k_B * this->TEMPERATURE;
    // calculate the partition function; Z(tau) == \Sum_{j=0}^{\Infinity} \exp{-E / \tau}
    for (unsigned int i = 0; i < this->states; i++) {
        this->P[i] = exp(-E.energy(i) / this->TAU);
        this->PARTITION += this->P[i];
    }
    // divide by Z to get the probability for each state
    for (unsigned int i = 0; i < this->states; i++) {
        this->P[i] /= this->PARTITION;
    }
}

/*
 * dynamically allocate the array and initialize everything
 * @param i             the size of the probability array / the number of states
 */
template <typename Num>
void PartitionFunctionSample<Num>::initialize(unsigned short int i) {
    this->P = new Num[i];
    this->states = i;
    for (unsigned int j = 0; j < this->states; j++) {
        this->P[j] = 0.0;
    }
}

/* class SystemManager */

template <typename Num>
SystemManager<Num>::~SystemManager(void) {
    delete [] this->sample;
}

template <typename Num>
SystemManager<Num>::SystemManager(void) {
    this->sample = nullptr;
}

template <typename Num>
void SystemManager<Num>::initialize(const std::string filename) {
    this->params.acquire(filename);
    this->sample = new PartitionFunctionSample<Num>[this->params.n_samp()];
    for (unsigned short i = 0; i < this->params.n_samp(); i++) {
        this->sample[i].initialize(this->params.states());
    }
}

template <typename Num>
bool SystemManager<Num>::save_to_disk(const std::string filename) {
    bool success = true;
    std::ofstream file(filename.c_str(), std::ofstream::out);
    if (file.is_open()) {
        // output the heading
        file << std::setprecision(16) << "All energies are in eV\n\nT (K),tau,Z(tau)";
        for (unsigned int i = 0; i < this->params.states(); i++) {
            file << ",P_" << i+1 << "(tau)";
        }
        file << '\n'; // start on the next row

        // output the data for each sample
        for (unsigned int i = 0; i < this->params.n_samp(); i++) {
            file << sample[i].T()   << ',' // temp
                 << sample[i].tau() << ',' // fundamental temp / thermal energy
                 << sample[i].Z();         // partition function
            for (unsigned int j = 0; j < this->params.states(); j++) {
                file << ',' << sample[i].P_i(j); // output the probabilities
            }
            file << '\n'; // next row
        }
    }
    else {
        success = false;
    }

    return success;
}


/* class progressBar */

template <typename Num>
progressBar<Num>::~progressBar(void) {
    this->full = this->current = 0;
    this->stream = nullptr;
}

template <typename Num>
progressBar<Num>::progressBar(void) {
    this->full = this->current = 0;
    this->stream = nullptr;
    this->width = 80;
}

template <typename Num>
progressBar<Num>::progressBar(unsigned int w) {
    this->full = this->current = 0;
    this->stream = nullptr;
    this->width = w;
}

/*
 * set up a text-based progress bar in the console output with a specified width
 * @param total         the full size of the index
 */
template <typename Num>
void progressBar<Num>::initialize(std::ostream& output, const Num total) {
    // set the 100% mark
    this->full = total;
    // save the pointer to the stream
    this->stream = &output;
    // print out the initial, empty indicator
    *(this->stream) << '[';
    for (unsigned int i = 0; i < (this->width); i++) {
        *(this->stream) << ' ';
    }
    // cap off the progress bar and return to the beginning of the line
    *(this->stream) << "]\r[|";
    // flush the buffer to make sure everything prints
    this->stream->flush();
}

/*
 * increments the progress bar by the fraction of the task that has been completed
 * since the last call
 * @param now           the current value of the index
 */
template <typename Num>
void progressBar<Num>::increment(const Num now) {
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

template <typename Num>
void progressBar<Num>::end(void) {
    *(this->stream) << '\n';
    this->stream->flush();
    this->current = static_cast<Num>(0);
    this->full = static_cast<Num>(0);
}

#endif
