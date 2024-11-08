#pragma once

#include "random.h"
#include "funzioni.h"
#include "wavesfunction.h"
#include <sstream>
#include <vector>

using namespace std;
using namespace arma;


const double _delta = 1.5; //parameter giving the jump to new position

//a_0 = 0.0529


// Abstract base class for the Metropolis algorithm
class System {

protected:
    
    double unif_delta = _delta;
    double pos_old = 0.;
    double _temp = 2.;
    double _energy_old = 0.753333333;
    double _err_old = 0.005;

    
    

public:
    Random RND;
    unique_ptr<Waves_function> WVT = unique_ptr<Wave_Variational_test>(new Wave_Variational_test()); //Pointer to wave function class
    unique_ptr<Waves_function> HT = unique_ptr<Energy_Variational_test>(new Energy_Variational_test()); //Pointer to Energy variational measure class

    // Constructor using initialization list to instantiate wave functions
    System();
    
    // Pure virtual methods to be implemented by derived classes
    

    double get_pos_old();
    void set_pos_new(double x);

    double get_delta();
    void set_delta(double x);

    double get_temp();
    void set_temp(double x);

    double get_energy_old();
    void set_energy_new(double x);

    double get_err_old();
    void set_err_new(double x);

    //metropolis method
    
    bool metropolis_method(Waves_function* WF); // Main Metropolis algorithm
    void initialize_acceptance(Waves_function* WF, double coordinates); // Initialize acceptance to target 50%

    void print_coordinate_accepted(double pos); //print in file.dat the coordinate accepted

    //simulated annealing method:

    void simulated_annealing (double temp, Waves_function* WF , Waves_function* Energy );

    bool Boltzmann(double temp, double Energy_old , double Energy_new);
    
    void cooling();

    void print_on_file (const string &filename,double content);
    void print_on_file (const string &filename,double mean, double err);
    void initialize_print();

    // Virtual destructor to properly delete dynamic wave function objects
     ~System() {}
};

// Derived class for Metropolis algorithm using Gaussian distribution







