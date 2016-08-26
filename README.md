# Partition Function Calculator
## Among other things

This program calculates the value of the partition function at a user-defined range of temperature points as well as the probability that a particle will be in each state defined in the partition function at each point. There are also tetration, factorial, and exponential functions accessible to the user, which are implemented by definition (tetration, factorial) or as a Taylor-Maclaurin series (e^x).

In order to compile this program, the Boost.Multiprecision and GMP libaries must be installed and linked (for g++, it's -lgmp and -lgmpxx) using the C++11 standard (-std=c++11 for g++).
