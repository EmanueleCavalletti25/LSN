#include <string>
#include "random.h"         
#include "funzioni.h"       
#include "wavesfunction.h"  
#include "system.h"    
#include <armadillo> 
#include <vector>

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



// Create unique pointers to the wave function and energy function tests
unique_ptr<Waves_function> WVT = unique_ptr<Wave_Variational_test>(new Wave_Variational_test());  // Pointer to wave function test
unique_ptr<Waves_function> HT = unique_ptr<Energy_Variational_test>(new Energy_Variational_test());  // Pointer to energy function test

WVT->set_mu(0.784);
WVT->set_sigma(0.618);


   

System sys;


sys.initialize_print();  // Initialize system

// Output streams for energy and positions
ofstream out_e, out_pos;
vector<double> energy(N_block, 0.0);  // Store energy per block
vector<double> energy_2(N_block, 0.0);  // Store squared energy per block

// Get mu and sigma for file naming
ostringstream oss_mu, oss_s;
oss_mu << fixed << std::setprecision(3) << WVT->get_mu();
oss_s << fixed << std::setprecision(3) << WVT->get_sigma();

// Initialize Metropolis algorithm and open output files
sys.initialize_acceptance(WVT.get(), pos_x_in);
//sys.set_delta(3.8+0.4*3.);
out_e.open("./OUTPUT/Energy_ave_" + oss_mu.str() + "_" + oss_s.str() + "_.dat");  // Energy output
out_pos.open("./OUTPUT/WaveFunction_" + oss_mu.str() + "_" + oss_s.str() + "_.dat", ios::app);  // Position output

// Equilibration loop (100,000 steps)
//for (int i = 0; i < 100000; i++) {
//    sys.metropolis_method(WVT.get());
//}



// Main loop over blocks

for (int j = 0; j < N_block; j++) {
    for (int i = 0; i < N_step; i++) {
        sys.metropolis_method(WVT.get());  // Execute Metropolis algorithm
        out_pos << sys.get_pos_old() << endl;  // Write position to output file
        
        energy[j] += HT->measure(sys.get_pos_old());  // Measure energy
        //cout << HT->measure(sys.get_pos_old()) << endl;
    }
    energy[j] /= (double)N_step;  // Average energy per block
    energy_2[j] = pow(energy[j], 2);  // Square the energy
}



// Perform block averaging and write results to energy file
Bloch_mean_method(energy, energy_2, N_block, out_e);

out_pos.close();  // Close position output file
return 0;
}
