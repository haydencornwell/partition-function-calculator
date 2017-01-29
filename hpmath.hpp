/*
 * High-precision mathematical functions
 */

#ifndef HPMATH_HPP
    #define HPMATH_HPP

#include <cmath>

//////////////////////////////
///// Function Templates /////
//////////////////////////////

/*
 * calculate a factorial of a nonnegative integer
 * @param n             the nonnegative integer
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

/*
 * calculate the natural logarithm of a positive real number
 * @param x             the number to calculate ln(x)
 * @return              ln(x)
 */
template <typename Numerical>
Numerical ln(const Numerical x) {
    Numerical y_n = x;
    Numerical y_np1 = 0.0;
    if (x > 0) {
        while (abs(y_np1 - y_n) > 1e-300) {
            y_n = y_np1;
            y_np1 = y_n + 2 * (x - exp(y_n)) / (x + exp(y_n));
        }
    }
    return y_np1;
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

#endif
