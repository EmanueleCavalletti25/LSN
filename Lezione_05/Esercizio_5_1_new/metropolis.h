#pragma once

#include "random.h"
#include "wavesfunction.h"
#include <sstream>

using namespace std;
using namespace arma;

const double a_0 = 0.0529 ; // Bohr radius constant

//a_0 = 0.0529


// Abstract base class for the Metropolis algorithm
class metropolis {

protected:
    double unif_sigma = a_0; // Sigma for uniform distribution
    double gauss_sigma = a_0; // Sigma for Gaussian distribution
    vec coord_unif_old = vec(3, fill::zeros); // Store old coordinates for uniform distribution
    vec coord_gauss_old = vec(3, fill::zeros); // Store old coordinates for Gaussian distribution
    

public:
    Random RND;
    unique_ptr<Waves_function> WF100 = unique_ptr<Waves_H_100>(new Waves_H_100()); // Pointer to wave function for 1s orbital
    unique_ptr<Waves_function> WF210 = unique_ptr<Waves_H_210>(new Waves_H_210()); // Pointer to wave function for 2p orbital
    
    
    // Constructor using initialization list to instantiate wave functions
    metropolis(){};
    
    // Pure virtual methods to be implemented by derived classes
    virtual double get_sigma() = 0;
    virtual void set_sigma(double x) = 0;
    virtual vec get_coord_old() = 0;
    virtual void set_coord_new(vec coord_new) = 0;
    virtual void metropolis_method(Waves_function* WF, vec coordinates) = 0; // Main Metropolis algorithm
    virtual void initialize_acceptance(Waves_function* WF, vec coordinates,int wf) = 0; // Initialize acceptance to target 50%
    virtual void print_coordinate_accepted(int wf,double pos) = 0; //print in file.dat the coordinate accepted
    

    // Virtual destructor to properly delete dynamic wave function objects
    virtual ~metropolis() {}
};

// Derived class for Metropolis algorithm using uniform distribution
class metropolis_unif : public metropolis {

public:
    metropolis_unif(); // Constructor

    // Implementations for uniform distribution-specific methods
    double get_sigma() override; // Get sigma for uniform distribution
    void set_sigma(double x) override; // Set sigma for uniform distribution
    vec get_coord_old() override ; // Get old coordinates for uniform distribution
    void set_coord_new(vec coord_new) override ; // Set new coordinates for uniform distribution

    void metropolis_method(Waves_function* WF, vec coordinates) override; // Implement Metropolis algorithm
    void initialize_acceptance(Waves_function* WF, vec coordinates,int wf) override; // Adjust parameters to achieve 50% acceptance
    void print_coordinate_accepted(int wf,double pos) override; // print a file.dat the coordinated accepted
};

// Derived class for Metropolis algorithm using Gaussian distribution
class metropolis_gauss : public metropolis {

public:
    metropolis_gauss(); // Constructor

    // Implementations for Gaussian distribution-specific methods
    double get_sigma() override; // Get sigma for Gaussian distribution
    void set_sigma(double x) override; // Set sigma for Gaussian distribution
    vec get_coord_old() override; // Get old coordinates for Gaussian distribution
    void set_coord_new(vec coord_new) override ; // Set new coordinates for Gaussian distribution

    void metropolis_method(Waves_function* WF, vec coordinates) override; // Implement Metropolis algorithm
    void initialize_acceptance(Waves_function* WF, vec coordinates,int wf) override; // Adjust parameters to achieve 50% acceptance
    void print_coordinate_accepted(int wf,double pos) override; // print a file.dat the coordinated accepted
};






