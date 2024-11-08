#include <string>
#include "random.h"         
#include "funzioni.h"       
#include "wavesfunction.h"  
#include "system.h"    
#include <armadillo> 
#include <vector>
#include <fstream>

using namespace std;
using namespace arma;       // Using Armadillo library for vectors, matrices
//int argc, char *argv[]
int main (int argc, char *argv[] ){

  if (argc < 2){
      cout << "Usage, select: " << argv [0] << "initial_position " << endl;
      return -1; 
   }

int N_block = 100;
int N_step = 100000;
double pos_x_in = atof(argv[1]);
double en_ave_new;
double en_err_new;
bool accept = false;

// Create unique pointers to the wave function and energy function tests
//unique_ptr<Waves_function> WVT = unique_ptr<Wave_Variational_test>(new Wave_Variational_test());  // Pointer to wave function test
//unique_ptr<Waves_function> HT = unique_ptr<Energy_Variational_test>(new Energy_Variational_test());  // Pointer to energy function test

System sys;


sys.initialize_print();  // Initialize system

// Output streams for energy and positions
ofstream out_e;
vector<double> energy(N_block, 0.0);  // Store energy per block
vector<double> energy_2(N_block, 0.0);  // Store squared energy per block

sys.HT->set_mu(0.);
sys.HT->set_sigma(0.5);
double mu_best = sys.HT->get_mu();
double sigma_best = sys.HT->get_sigma();
int contatore = 1;
double en_t = 0.;
double mu_t = 0.;
double s_t = 0.;

// Begin the annealing schedule
while (sys.get_temp() > 0.01) {
    
    // Loop for 10 iterations to calculate energy
    while (contatore <= 3) {
        cout << "Cycle number: " << contatore << endl;     
        vector<double> energy(N_block, 0.0);  // Store energy per block
        vector<double> energy_2(N_block, 0.0);  // Store squared energy per block

        // Print current mu and sigma values to file
        sys.print_on_file("./OUTPUT/mu_sigma_proposti.dat", sys.HT->get_mu(), sys.HT->get_sigma());

        // Set delta for Metropolis algorithm
        sys.set_delta(sys.WVT->get_mu() + 3 * sys.WVT->get_sigma());

        // Main loop over blocks
        for (int j = 0; j < N_block; j++) {
            for (int i = 0; i < N_step; i++) {
                sys.metropolis_method(sys.WVT.get());  // Execute Metropolis algorithm
                // Ensure mu value has not changed unexpectedly
                if (sys.WVT->get_mu() != sys.HT->get_mu()) {
                    exit(-1);  // Exit if mu has changed
                }
                energy[j] += sys.HT->measure(sys.get_pos_old());  // Measure energy
            }
            energy[j] /= (double)N_step;  // Average energy per block
            energy_2[j] = pow(energy[j], 2);  // Square the energy
        }

        // Perform block averaging and write results to energy file
        Bloch_mean_method(energy, energy_2, N_block, en_ave_new, en_err_new);

        // Calculate acceptance using the Boltzmann distribution
        accept = sys.Boltzmann(sys.get_temp(), sys.get_energy_old(), en_ave_new);

        // If accepted, update energy and best parameters
        if (accept == true) {
            sys.set_energy_new(en_ave_new);
            sys.set_err_new(en_err_new);
            pos_x_in = 0.0;
            mu_best = sys.HT->get_mu();
            sigma_best = sys.HT->get_sigma();

            cout << "Accepted step!" << endl << "Temperature: " << sys.get_temp() << "  Energy: " << en_ave_new << " Error: " << en_err_new << endl;
        }
        // If not accepted, revert to the best parameters
        if (accept == false) {
            sys.HT->set_mu(mu_best);
            sys.HT->set_sigma(sigma_best);
            sys.set_energy_new(sys.get_energy_old());

            cout << "Rejected step!" << endl << "Temperature: " << sys.get_temp() << "  Energy: " << sys.get_energy_old() << " Error: " << sys.get_err_old() << endl;
        }
        
        // Update mu and sigma with a random adjustment
        sys.WVT->set_mu(abs(sys.WVT->get_mu() + sys.get_temp() * sys.RND.Rannyu(-1, 1) * 0.5));
        sys.WVT->set_sigma(abs(sys.WVT->get_sigma() + sys.get_temp() * sys.RND.Rannyu(-1, 1) * 0.25));
        contatore++;  // Increment the counter
        pos_x_in = 0;  // Reset position variable
    }
    
    // Print average energy and best parameters to output files
    sys.print_on_file("./OUTPUT/Energy_ave.dat", sys.get_energy_old(), sys.get_err_old());
    sys.print_on_file("./OUTPUT/mu_sigma.dat", mu_best, sigma_best);
    
    // Cool the system
    sys.cooling();
    contatore = 1;  // Reset counter for the next temperature stage

    // Check if the current energy is the best
    if (sys.get_energy_old() < en_t) {
        en_t = sys.get_energy_old();  // Update best energy
        mu_t = sys.WVT->get_mu();     // Update best mu
        s_t = sys.WVT->get_sigma();    // Update best sigma
    }
}

// Output the best measured values
ofstream out;
out.open("./OUTPUT/Best_measure.dat");
out << "Best Energy" << setw(12) << "Best MU" << setw(12) << "Best SIGMA" << endl;
out << en_t << setw(12) << mu_t << setw(12) << s_t << endl;

cout << "Best Energy" << setw(12) << "Best MU" << setw(12) << "Best SIGMA" << endl;
cout << en_t << setw(12) << mu_t << setw(12) << s_t << endl;

cout << "Annealing finished!!" << endl;



  // Close position output file
return 0;
}



//ciclo; proposta mu e sigma e aggiornare con quelli vecchi!

