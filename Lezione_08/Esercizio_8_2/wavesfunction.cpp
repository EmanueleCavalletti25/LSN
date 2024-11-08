#include "wavesfunction.h"

// Definition of global parameters for wave function
double _mu = 0.0;   // Parameter mu for wave function
double _sigma = 1.0; // Parameter sigma for wave function

// Initialization of static members in Waves_function class
double Waves_function::unif_mu = _mu;  // Initialize static mu
double Waves_function::unif_sigma = _sigma;  // Initialize static sigma

// Constructor for the base class Waves_function
Waves_function::Waves_function() {};

// Virtual destructor for the base class Waves_function
Waves_function::~Waves_function() {};

// Constructor for the derived class Wave_Variational_test
Wave_Variational_test::Wave_Variational_test() {};

// Implementation of the measure function for Wave_Variational_test
double Wave_Variational_test::measure(double x) {
    // Returns the sum of two Gaussian wave functions (one shifted by +mu and one by -mu)
    return exp(-1. * pow((x - unif_mu), 2) / (2 * pow(unif_sigma, 2))) + 
           exp(-1. * pow((x + unif_mu), 2) / (2 * pow(unif_sigma, 2)));
}

// Getter for sigma (shared among all derived classes)
double Wave_Variational_test::get_sigma() {
    return unif_sigma;  // Returns the shared static sigma value
}

// Setter for sigma (modifies the shared sigma for all derived classes)
void Wave_Variational_test::set_sigma(double x) {
    unif_sigma = x;  // Modifies the shared static sigma value
}

// Getter for mu (shared among all derived classes)
double Wave_Variational_test::get_mu() {
    return unif_mu;  // Returns the shared static mu value
}

// Setter for mu (modifies the shared mu for all derived classes)
void Wave_Variational_test::set_mu(double x) {
    unif_mu = x;  // Modifies the shared static mu value
}

// Implementation of the Energy_Variational_test derived class
Energy_Variational_test::Energy_Variational_test() {};

// Measure function for Energy_Variational_test
double Energy_Variational_test::measure(double x) {
    // Create a unique pointer to the base class Waves_function (specifically Wave_Variational_test)
   

    // Measure the wave function at position x
    double wave = exp(-1. * pow((x - unif_mu), 2) / (2 * pow(unif_sigma, 2))) + 
           exp(-1. * pow((x + unif_mu), 2) / (2 * pow(unif_sigma, 2)));


    // Define the potential function (for example, x^4 - (5/2)*x^2)
    double potential = pow(x, 4) - (5. / 2.) * pow(x, 2);

    // Calculate the second derivative of the wave function with respect to x
    double second_derivate = 
        (-1. / pow(unif_sigma, 2) + pow((x + unif_mu), 2) / pow(unif_sigma, 4)) * 
        exp(-1. * pow((x + unif_mu), 2) / (2. * pow(unif_sigma, 2))) + 
        (-1. / pow(unif_sigma, 2) + pow((x - unif_mu), 2) / pow(unif_sigma, 4)) * 
        exp(-1. * pow((x - unif_mu), 2) / (2. * pow(unif_sigma, 2)));

    // Return the energy expectation value (Laplacian term + potential)
    return (-0.5 * (second_derivate / wave) + potential);
}

// Getter for sigma in Energy_Variational_test (shared across all derived classes)
double Energy_Variational_test::get_sigma() {
    return unif_sigma;  // Returns the shared static sigma value
}

// Setter for sigma in Energy_Variational_test (modifies the shared static sigma)
void Energy_Variational_test::set_sigma(double x) {
    unif_sigma = x;  // Modifies the shared static sigma value
}

// Getter for mu in Energy_Variational_test (shared across all derived classes)
double Energy_Variational_test::get_mu() {
    return unif_mu;  // Returns the shared static mu value
}

// Setter for mu in Energy_Variational_test (modifies the shared static mu)
void Energy_Variational_test::set_mu(double x) {
    unif_mu = x;  // Modifies the shared static mu value
}
