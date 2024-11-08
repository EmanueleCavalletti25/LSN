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

// Default global parameters for the wave function
extern double _mu;  // Mu parameter for the wave function
extern double _sigma;  // Sigma parameter for the wave function

// Abstract base class for wave functions
class Waves_function {

    // Static members that will be shared across all derived classes
    protected:
    static double unif_mu;   // Mu shared across all derived classes
    static double unif_sigma;  // Sigma shared across all derived classes

public:
    Waves_function(); // Constructor

    // Pure virtual functions that need to be implemented by derived classes
    virtual double get_sigma() = 0;  // Get sigma
    virtual void set_sigma(double x) = 0;  // Set sigma

    virtual double get_mu() = 0;  // Get mu
    virtual void set_mu(double x) = 0;  // Set mu

    virtual double measure(double x) = 0;  // Pure virtual function for measurement
    virtual ~Waves_function();  // Virtual destructor
};

// Derived class for the variational wave function test
class Wave_Variational_test : public Waves_function {
public:
    Wave_Variational_test();  // Constructor

    double measure(double x) override;  // Implementation of the wave variational test

    double get_sigma() override;  // Get sigma
    void set_sigma(double x) override;  // Set sigma

    double get_mu() override;  // Get mu
    void set_mu(double x) override;  // Set mu
};

// Derived class for the energy variational test
class Energy_Variational_test : public Waves_function {
public:
    Energy_Variational_test();  // Constructor
    
    double measure(double x) override;  // Implementation of the energy LOC (Laplacian Operator Calculation)

    double get_sigma() override;  // Get sigma
    void set_sigma(double x) override;  // Set sigma

    double get_mu() override;  // Get mu
    void set_mu(double x) override;  // Set mu
};
