#ifndef HPMATH_HPP
    #define HPMATH_HPP

#include <iostream>
    using std::cin;
    using std::cout;
    using std::endl;
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
    ~partition_fxn_sample(void);
    partition_fxn_sample(const unsigned int);
    partition_fxn_sample(void);
    mpf_float_1000 tau; // fundamental temperature
    mpf_float_1000 Z; // partition function at tau
    mpf_float_1000* P; // owning pointer to state probability array
    unsigned short int n; // size of state probability array
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
};

///////////////////////////////////////
///// Member Function Definitions /////
///////////////////////////////////////

/* class partition_fxn_sample */

partition_fxn_sample::~partition_fxn_sample(void) {
    // clean up
    delete [] this->P;
    this->P = nullptr;
    this->n = 0;
}

partition_fxn_sample::partition_fxn_sample(void) {
    // empty
    this->tau = this->Z = 0.0;
    this->P = nullptr;
    this->n = 0;
}

partition_fxn_sample::partition_fxn_sample(const unsigned int numstates) {
    // allocate the probability state array
    this->P = new mpf_float_1000[numstates];
    this->n = numstates;
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

template <typename Num>
void progressBar<Num>::initialize(std::ostream& output, const Num total) {
    // set the 100% mark
    this->full = total;
    // save the pointer to the stream
    this->stream = &output;
    // print out the initial, empty indicator
    *(this->stream) << '[';
    for (unsigned int i = 0; i < (this->width-1); i++) {
        *(this->stream) << ' ';
    }
    // cap off the progress bar and return to the beginning of the line
    *(this->stream) << "]\r[|";
    // flush the buffer to make sure everything prints
    this->stream->flush();
}

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

//////////////////////////////
///// Function Templates /////
//////////////////////////////

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
 * A template for a tetration function which calculates a power tower
 * and repeated roots for negative hyperpowers
 * @ param base         the number to be tetrated
 * @ param hyperpower   the height of the power/root tower
 * @ return             the result of the tetration
 */
template <typename Numerical>
Numerical tetrate(const Numerical base, const int hyperpower) {
    Numerical result = base;
    // repeated exponentiation for positive hyperpowers
    if (hyperpower > 0) {
        for (int i = 1; i < hyperpower; i++) {
            result = pow(base, result);
        }
    }
    // repeated roots for negative powers
    else if (hyperpower < 0 && base > 0) {
        // decrement since the hyperpower is negative
        for (int i = 1; i < -hyperpower; i++) {
            result = pow(result, 1.0 / base);
        }
    }
    // return 1.00 for zero powers
    else {
        result = 1;
    }

    return result;
}

/*
 * Raises Euler's number to a power (doesn't stop until all digits are filled!)
 * @param x             the number to raise e to
 * @return              the result
 */
template <typename Numerical>
Numerical exp(const Numerical x) {
    Numerical result = 1; // 0th term == 1
    Numerical prev = 0;
    unsigned long long int i = 1;
    // e^x == \Sum_{n==0}^{Inf} \frac{x^n}{n!}
    while (result != prev && i <= std::numeric_limits<unsigned long long int>::max()) {
        prev = result;
        Numerical numerator = 1;
        // x^i on index i without using double pow(mpf_float_1000, mpz_int)
        for (unsigned long long int j = 1; j <= i; j++) {
            numerator *= x;
        }
        result += numerator / factorial(i);
        i++;
    }

    return result;
}

#endif