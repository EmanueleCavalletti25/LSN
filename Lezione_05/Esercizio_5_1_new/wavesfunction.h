#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <armadillo>
#include <cstdlib>
#include <vector>
#include <memory>
#include <iomanip>

using namespace std;
using namespace arma;



// Abstract base class for wave functions
class Waves_function {
public:
    Waves_function(); // Constructor

    virtual double measure(vec coordinates) = 0; // Pure virtual function for measurement
    virtual ~Waves_function(); // Virtual destructor
};

// Wave function for Hydrogen atom's 1s orbital (n=1, l=0, m=0)
class Waves_H_100 : public Waves_function {
public:
    Waves_H_100(); // Constructor

    double measure(vec coordinates) override; // Implementation of 1s wave function measurement
};

// Wave function for Hydrogen atom's 2p orbital (n=2, l=1, m=0)
class Waves_H_210 : public Waves_function {
public:
    Waves_H_210(); // Constructor

    double measure(vec coordinates) override; // Implementation of 2p wave function measurement
};