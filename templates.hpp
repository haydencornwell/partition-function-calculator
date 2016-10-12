#ifndef TEMPLATES_HPP
    #define TEMPLATES_HPP

#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include <fstream>
#include <cmath>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/cpp_int.hpp>
    using namespace boost::multiprecision;

/////////////////////
///// Operators /////
/////////////////////

mpfr_float_1000 operator"" _mpr1k(long double x) {return mpfr_float_1000(x);}
mpfr_float_1000 operator"" _mpr1k(unsigned long long int x) {return mpfr_float_1000(x);}

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
