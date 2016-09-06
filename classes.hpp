#ifndef CLASSES_HPP
    #define CLASSES_HPP

/////////////////
///// Enums /////
/////////////////

enum MenuPick {
    tetration,
    factorials,
    exponentiate,
    partition,
    log_base_n,
    quit
};

///////////////////
///// Objects /////
///////////////////

class partition_fxn_sample {
    const mpfr_float_1000 k_B = 8.6173324e-5; // (eV/K) Boltzmann constant
    unsigned short int n; // size of state probability array
    mpfr_float_1000 tau; // fundamental temperature
    mpfr_float_1000 Z; // partition function at tau
    mpfr_float_1000* P; // owning pointer to state probability array
    mpfr_float_1000 T; // temperature in K
  public:
    ~partition_fxn_sample(void);
    partition_fxn_sample(const unsigned int);
    partition_fxn_sample(void);
    // accessors
    mpfr_float_1000 get_tau(void) {return this->tau;}
    mpfr_float_1000 get_Z(void) {return this->Z;}
    mpfr_float_1000 get_P_i(unsigned short int i) {return (i < n ? this->P[i] : 0);}
    mpfr_float_1000 get_T(void) {return this->T;}
    // calculation
    template <typename Numerical>
    void calculate(Numerical, Numerical*);
    // initialization in case the default constructor was used
    void initialize(unsigned short int);
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
    this->tau = this->Z = 0.0;
    // allocate the probability state array
    this->P = new mpfr_float_1000[numstates];
    this->n = numstates;
    for (unsigned int i = 0; i < this->n; i++) {
        this->P[i] = 0.0;
    }
}

/*
 * calculate the values at the given temperature and state energies
 * @param temp          the current temperature in Kelvins
 * @param E             a pointer to an array of state energy values in eV
 */
template <typename Numerical>
void partition_fxn_sample::calculate(Numerical temp, Numerical* E) {
    this->T = temp;
    this->tau = this->k_B * this->T;
    // calculate the partition function; Z(tau) == \Sum_{j=0}^{\Infinity} \exp{-E / \tau}
    for (unsigned int i = 0; i < this->n; i++) {
        this->P[i] = exp(-E[i] / this->tau);
        this->Z += this->P[i];
    }
    // divide by Z to get the probability for each state
    for (unsigned int i = 0; i < this->n; i++) {
        this->P[i] /= this->Z;
    }
}

/*
 * dynamically allocate the array and initialize everything
 * @param i             the size of the array
 */
void partition_fxn_sample::initialize(unsigned short int i) {
    this->P = new mpfr_float_1000[i];
    this->n = i;
    for (unsigned int j = 0; j < this->n; j++) {
        this->P[j] = 0.0;
    }
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

#endif
