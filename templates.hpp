/*
 * Assorted free-floating functions and templates
 */

#ifndef TEMPLATES_HPP
    #define TEMPLATES_HPP

#include <iostream>
#include <limits>
#include <string>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/mpfr.hpp>
    using namespace boost::multiprecision;

/////////////////////
///// Operators /////
/////////////////////

mpfr_float_1000 inline operator"" _mpr1k(long double x) {return mpfr_float_1000(x);}
mpfr_float_1000 inline operator"" _mpr1k(unsigned long long int x) {return mpfr_float_1000(x);}

/////////////////////
///// Constants /////
/////////////////////
/**
 * mathematical and physical constants, as mpfr_float_1000's
 */
namespace Constants {
    //! (eV/K) Boltzmann constant
    const mpfr_float_1000 k_B = 8.6173324e-5_mpr1k;
    //! (per mol) Avogadro constant
    const mpfr_float_1000 N_A = 6.022140857e23_mpr1k;
    //! (V) absolute potential of an electron at rest in a vacuum vs SHE
    const mpfr_float_1000 V_abs = 4.44_mpr1k;
    //! (unitless) pi, to 1000 digits
    const mpfr_float_1000 mpfr_pi = boost::math::constants::pi<mpfr_float_1000>();
    //! (F / m) electric permittivity
    const mpfr_float_1000 permittivity = 8.854187817e-12_mpr1k;
    //! (N * m^2 / C^2) Coulomb's constant
    const mpfr_float_1000 k_e = 1_mpr1k / (4_mpr1k * mpfr_pi * permittivity);
}

///////////////////////////////////
///// Free-floating functions /////
///////////////////////////////////

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
 * get checked input until the input is successful
 * @param in            the stream to read from
 * @param out           the stream to write error messages to
 * @param val           the variable to save to
 * @param err_msg       the error message to display if the input is rejected
 */
template <typename T>
void getterLoop(std::istream& in, std::ostream& out, T& val, std::string err_msg) {
    bool good;
    do {
        good = getInput<T>(in, val);
        if (!good) {
            out << err_msg;
        }
    } while (!good);
}

/*
 * get checked input until the input is successful and it falls within a given range
 * @param in            the stream to read from
 * @param out           the stream to write error messages to
 * @param val           the variable to save to
 * @param min           the minimum acceptable value (inclusive)
 * @param max           the maximum acceptable value (inclusive)
 * @param err_msg       the error message to display if the input is rejected
 */
template <typename T>
void rangedGetterLoop(std::istream& in, std::ostream& out, T& val, const T min, const T max, std::string err_msg) {
    bool good;
    do {
        good = getRangedInput<T>(in, val, min, max);
        if (!good) {
            out << err_msg;
        }
    } while (!good);
}

#endif
