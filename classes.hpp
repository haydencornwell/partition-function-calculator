/**
 * Class definition file
 */

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

///////////////////
///// Objects /////
///////////////////

/**
 * encapsulates the classes describing a thermodynamic system
 */
namespace Thermodynamics{
    /**
     * acquires and contains the parameters for a system
     */
    template <typename Num>
    class SystemParameters {
        //! size of the energy array
        unsigned short int n;
        //! (eV) owning pointer to energy array
        Num* E;
        //! (eV) owning pointer to array of total chemical potentials
        Num* TOTAL_POTENTIAL;
        //! (K) temperature
        Num TEMPERATURE;
      public:
        //! the name of the file to save to
        std::string filename;
        ~SystemParameters(void);
        SystemParameters(void);
        // accessors
        //! return the number of states in the system
        unsigned short int states(void) const {return this->n;}
        //! return the temperature of the system
        Num T(void) const {return this->TEMPERATURE;}
        void set_T(const Num temp) {this->TEMPERATURE = temp;}
        void set_mu(const unsigned short i, const Num potential) {this->TOTAL_POTENTIAL[i] = potential;}
        //! return the total chemical potential of a state
        Num mu(unsigned short int i) const {return (i < this->n ? this->TOTAL_POTENTIAL[i] : static_cast<Num>(0));}
        //! return the energy of a state
        Num energy(unsigned short int i) const {return (i < this->n ? this->E[i] : static_cast<Num>(0));}
        // other functions
        SystemParameters& acquire(const std::string);
    };

    /**
     * contains sample data and has the ability to calculate its value from a SystemParameters object
     */
    template <typename Num>
    class PartitionFunctionSample {
        //! size of state probability array
        unsigned short int states;
        //! (eV) fundamental temperature
        Num TAU;
        //! partition function at tau
        Num PARTITION;
        //! owning pointer to state probability array
        Num* P;
        //! (eV) owning pointer to array of total chemical potentials
        Num* TOTAL_POTENTIAL;
        //! (K) temperature
        Num TEMPERATURE;
      public:
        ~PartitionFunctionSample(void);
        PartitionFunctionSample(const unsigned int);
        PartitionFunctionSample(void);
        // accessors
        //! return the fundamental temperature (thermal energy) of the system
        Num tau(void) const {return this->TAU;}
        //! return the value of the partition function
        Num Z(void) const {return this->PARTITION;}
        //! return the probability of the given state
        Num P_i(unsigned short int i) const {return (i < this->states ? this->P[i] : static_cast<Num>(0));}
        //! return the chemical potential of the given state
        Num mu_i(unsigned short int i) const {return (i < this->states ? this->TOTAL_POTENTIAL[i] : static_cast<Num>(0));}
        //! return the temperature of the system
        Num T(void) const {return this->TEMPERATURE;}
        // calculation
        void calculate(SystemParameters<Num>&);
        // initialization in case the default constructor was used
        void initialize(unsigned short int);
    };

    /**
     * manages sample data and loading/saving for a single system
     */
    template <typename Num>
    class SystemManager {
        unsigned short int number_of_samples = 0;
      public:
        //! thermodynamic system parameters
        SystemParameters<Num> params;
        //! owning pointer to sample array
        PartitionFunctionSample<Num>* sample;
        ~SystemManager(void);
        SystemManager(void);
        unsigned short int n_samp(void) {return number_of_samples;}
        bool save_to_disk(std::string);
        void initialize(const unsigned short int);
    };
}

/**
 * draws a simple progress bar in the console with text
 */
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

////////////////////////////
/* class SystemParameters */

/**
 * destructor
 */
template <typename Num>
Thermodynamics::SystemParameters<Num>::~SystemParameters(void) {
    delete [] this->E;
    delete [] this->TOTAL_POTENTIAL;
    this->E = nullptr;
    this->TOTAL_POTENTIAL = nullptr;
    this->n = 0;
}

/**
 * default constructor
 */
template <typename Num>
Thermodynamics::SystemParameters<Num>::SystemParameters(void) {
    this->E = nullptr;
    this->TOTAL_POTENTIAL = nullptr;
    this->n = 0;
}

/**
 * acquire system information from a file or user input
 * @param cfg_name      the name of the config file
 * @return              itself, by reference
 */
template <typename Num>
Thermodynamics::SystemParameters<Num>& Thermodynamics::SystemParameters<Num>::acquire(const std::string cfg_name) {
    std::ifstream config(cfg_name);
    std::string use_cfg_response;
    // try to read from a config file first
    if (config.is_open()) {
        std::cout << "\nConfiguration file found; use data? (y/n) ";
        std::getline(std::cin, use_cfg_response);
        if (static_cast<char>(tolower(use_cfg_response[0])) == 'y') {
            std::getline(config, this->filename, '\n');
            // handle Windows/DOS line endings when using std::getline
            if (this->filename.back() == '\r') {
                this->filename.pop_back(); // delete the last character, which is a carriage return
            }
            config >> this->n;
            if (this->n != 0) {
                this->E = new Num[this->n];
                this->TOTAL_POTENTIAL = new Num[this->n];
                for (unsigned short int i = 0; i < this->n; i++) {
                    config >> this->E[i] >> this->TOTAL_POTENTIAL[i];
                }
            }
            config.close();
        }
        else {
            config.close();
        }

        if (!config.good() || this->n == 0) {
            std::cout << "Error reading configuration file.\n\n";
            config.close();
        }
    }

    if (static_cast<char>(tolower(use_cfg_response[0])) == 'n' || !config.good()) {
        // file name
        std::cout << "\nEnter a filename to save the results (CSV format, will be overwritten): ";
        std::getline(std::cin, this->filename);

        std::cout << "How many states does the partition function have? ";
        rangedGetterLoop(std::cin, std::cout, this->n, static_cast<unsigned short>(0),
                         std::numeric_limits<unsigned short int>::max(),
                         "Please enter a positive integer: ");

        this->E = new Num[this->n];
        this->TOTAL_POTENTIAL = new Num[this->n];

        // get the energies
        for (unsigned short int i = 0; i < this->n; i++) {
            std::cout << "Enter the energy of the " << i+1
                      << ((i+1 % 10 == 1 && i+1 % 100 != 11) ? "st" :
                      ((i+1 % 10 == 2 && i+1 % 100 != 12) ? "nd" :
                      ((i+1 % 10 == 3 && i+1 % 100 != 13) ? "rd" : "th")))
                      << " state in eV: ";
            getterLoop(std::cin, std::cout, this->E[i], "Please enter a numerical value: ");
        }
    }

    // returns itself so this can be called inside of another function that uses this type.
    return *this;
}

///////////////////////////////////
/* class PartitionFunctionSample */

/**
 * destructor
 */
template <typename Num>
Thermodynamics::PartitionFunctionSample<Num>::~PartitionFunctionSample(void) {
    delete [] this->P;
    delete [] this->TOTAL_POTENTIAL;
    this->P = nullptr;
    this->TOTAL_POTENTIAL = nullptr;
    this->states = 0;
}

/**
 * default constructor
 */
template <typename Num>
Thermodynamics::PartitionFunctionSample<Num>::PartitionFunctionSample(void) {
    this->TAU = this->PARTITION = 0.0;
    this->P = nullptr;
    this->TOTAL_POTENTIAL = nullptr;
    this->states = 0;
}

/**
 * constructor with dynamic allocation
 */
template <typename Num>
Thermodynamics::PartitionFunctionSample<Num>::PartitionFunctionSample(const unsigned int numstates) {
    this->TAU = this->PARTITION = 0.0;
    this->P = new Num[numstates];
    this->TOTAL_POTENTIAL = new Num[numstates];
    this->states = numstates;
    for (unsigned int i = 0; i < this->states; i++) {
        this->P[i] = this->TOTAL_POTENTIAL = 0.0;
    }
}

/**
 * calculate the values at the given temperature and state energies
 * @param E             a pointer to system parameters
 */
template <typename Num>
void Thermodynamics::PartitionFunctionSample<Num>::calculate(SystemParameters<Num>& params) {
    this->TEMPERATURE = params.T();
    this->TAU = Constants::k_B * this->TEMPERATURE;
    // calculate the partition function; Z(tau) == \Sum_{j=0}^{\Infinity} \exp{(\mu - E) / \tau}
    for (unsigned int i = 0; i < this->states; i++) {
        this->P[i] = exp((params.mu(i) - params.energy(i)) / this->TAU);
        this->PARTITION += this->P[i];
        // this is just for bookkeeping purposes
        this->TOTAL_POTENTIAL[i] = params.mu(i);
    }
    // divide by Z to get the probability for each state
    for (unsigned int i = 0; i < this->states; i++) {
        this->P[i] /= this->PARTITION;
    }
}

/**
 * dynamically allocate the array and initialize everything
 * @param i             the size of the probability array / the number of states
 */
template <typename Num>
void Thermodynamics::PartitionFunctionSample<Num>::initialize(const unsigned short int i) {
    this->P = new Num[i];
    this->TOTAL_POTENTIAL = new Num[i];
    this->states = i;
    for (unsigned int j = 0; j < this->states; j++) {
        this->P[j] = this->TOTAL_POTENTIAL[j] = 0.0;
    }
}

/////////////////////////
/* class SystemManager */

/**
 * destructor
 */
template <typename Num>
Thermodynamics::SystemManager<Num>::~SystemManager(void) {
    delete [] this->sample;
}

/**
 * default constructor
 */
template <typename Num>
Thermodynamics::SystemManager<Num>::SystemManager(void) {
    this->sample = nullptr;
}

/**
 * allocate the sample object array and initalize each element
 * @param n_samp        the number of samples
 */
template <typename Num>
void Thermodynamics::SystemManager<Num>::initialize(const unsigned short int n_samp) {
    this->sample = new PartitionFunctionSample<Num>[n_samp];
    this->number_of_samples = n_samp;
    for (unsigned short i = 0; i < n_samp; i++) {
        this->sample[i].initialize(this->params.states());
    }
}

/**
 * save the results to disk
 * @param filename      The name of the save file
 * @return              whether or not the save was successful
 */
template <typename Num>
bool Thermodynamics::SystemManager<Num>::save_to_disk(const std::string filename) {
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
        for (unsigned int i = 0; i < this->n_samp(); i++) {
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

///////////////////////
/* class progressBar */

/**
 * destructor
 */
template <typename Num>
progressBar<Num>::~progressBar(void) {
    this->full = this->current = 0;
    this->stream = nullptr;
}

/**
 * default constructor
 */
template <typename Num>
progressBar<Num>::progressBar(void) {
    this->full = this->current = 0;
    this->stream = nullptr;
    this->width = 80;
}

/**
 * constructor with a custom width
 * @param w             the custom progress bar width, in characters
 */
template <typename Num>
progressBar<Num>::progressBar(unsigned int w) {
    this->full = this->current = 0;
    this->stream = nullptr;
    this->width = w;
}

/**
 * set up a text-based progress bar in the console output with a specified width
 * @param output        the output stream to write to
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

/**
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

/**
 * end the progress bar, writes a newline, and resets the progress bar's internal state
 */
template <typename Num>
void progressBar<Num>::end(void) {
    *(this->stream) << '\n';
    this->stream->flush();
    this->current = static_cast<Num>(0);
    this->full = static_cast<Num>(0);
}

#endif
