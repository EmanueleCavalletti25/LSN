#include <string>
#include "random.h"         
#include "funzioni.h"       
#include "wavesfunction.h"  
#include "metropolis.h"     


using namespace std;
using namespace arma;       // Using Armadillo library for vectors, matrices

int main ( int argc, char *argv[]){
  // Check command-line arguments
  if (argc < 3){
      cout << "Usage, select: " << argv [0] << "<method: gaussian/uniform> wave (n,l,m): <100 or 210>"  << endl;
      return -1;
  }

  string method = argv[1];  // Choose method: "gaussian" or "uniform"
  int wf = stoi(argv[2]);   // Select quantum number: 100 or 210

  // Verify if the quantum number is correct
  if (wf != 100 && wf != 210 ){
    cerr << "write the quantic number (n,l,m): 100 or 210" << endl;
    return -1;
  }

  // Initial coordinates (X=50, Y=Z=0)
  vec in_coordinates(3,fill::zeros); 
  in_coordinates[0]= 0.;
  double N_block = 500;      // Number of blocks
  double step_block = 1000;  // Number of steps per block

  // Create strings to save output files
  ostringstream oss_x, oss_y, oss_z;
  oss_x << setprecision(1) << fixed << in_coordinates[0];
  oss_y << setprecision(1) << fixed << in_coordinates[1];
  oss_z << setprecision(1) << fixed << in_coordinates[2];

  // Loop to vary the number of steps per block
  for(int l=0; l<12; l++){  //remove this loop to obtain only one cicle, and select the step. Remove also the curly bracket at line 150
    ostringstream oss_step;
    vec r(step_block, fill::zeros);    // Vector for average radius
    vec r_2(N_block, fill::zeros);     // Vector for square of average radius
    vec r_blk(N_block, fill::zeros);   // Block for mean calculations
    step_block = 1000 + l*1000;        // Increment steps per block
    oss_step << setprecision(0) << fixed << step_block;

    cout << "number of passes for block: " << step_block << endl;

    // Check if the method is correct
    if (method != "gaussian" && method != "uniform" ){
      cerr << "write the method: gaussian or uniform" << endl;
      return -1;
    } 

    // Gaussian method
    else if (method == "gaussian"){
      // Initialize the metropolis object for Gaussian method
      unique_ptr<metropolis> gauss_metropolis = unique_ptr<metropolis_gauss>(new metropolis_gauss);
      gauss_metropolis->set_coord_new(in_coordinates);

      // If quantum number is 100
      if( wf == 100){
        ofstream out_g_100;
        gauss_metropolis->initialize_acceptance(gauss_metropolis->WF100.get(), in_coordinates, wf);  // Initialize Metropolis
        out_g_100.open("./OUTPUT/raggio_gauss_w100_pos_"+oss_x.str()+"_"+oss_y.str()+"_"+oss_z.str()+"_"+oss_step.str()+"_step.dat");
        for (int i = 0; i < 100000; i++){
          gauss_metropolis->metropolis_method(gauss_metropolis->WF100.get(), in_coordinates);
        }
        // Loop over blocks
        for (int j = 0; j < N_block; j++){
          for (int i = 0; i < step_block; i++) {
            gauss_metropolis->metropolis_method(gauss_metropolis->WF100.get(), in_coordinates);  // Execute Metropolis
            vec radius = gauss_metropolis->get_coord_old();  // Get current coordinates
            r[j] += radius[0];  // Add radius
          }
          r[j] /= step_block;       // Calculate mean radius per block
          r_2[j] = pow(r[j],2);     // Calculate squared mean radius
        }
        // Perform block averaging and write results
        Bloch_mean_method(r, r_2, N_block, out_g_100);
      }

      // If quantum number is 210
      else if( wf == 210){
        ofstream out_g_210;
        gauss_metropolis->set_coord_new(in_coordinates);
        gauss_metropolis->initialize_acceptance(gauss_metropolis->WF210.get(), in_coordinates, wf);  // Initialize Metropolis for 210
        out_g_210.open("./OUTPUT/raggio_gauss_w210_pos_"+oss_x.str()+"_"+oss_y.str()+"_"+oss_z.str()+"_"+oss_step.str()+"_step.dat");

        for (int i = 0; i < 100000; i++){
          gauss_metropolis->metropolis_method(gauss_metropolis->WF100.get(), in_coordinates);
        }
        // Loop over blocks
        for (int j = 0; j < N_block; j++){
          for (int i = 0; i < step_block; i++) {
            gauss_metropolis->metropolis_method(gauss_metropolis->WF210.get(), in_coordinates);  // Execute Metropolis
            vec radius = gauss_metropolis->get_coord_old();  // Get current coordinates
            r[j] += radius[0];  // Add radius
          }
          r[j] /= step_block;       // Calculate mean radius per block
          r_2[j] = pow(r[j],2);     // Calculate squared mean radius
        }
        // Perform block averaging and write results
        Bloch_mean_method(r, r_2, N_block, out_g_210);
      }
    }

    // Uniform method
    else if (method == "uniform"){
      // Initialize the metropolis object for uniform method
      unique_ptr<metropolis> unif_metropolis = unique_ptr<metropolis_unif>(new metropolis_unif);
      unif_metropolis->set_coord_new(in_coordinates);

      // If quantum number is 100
      if( wf == 100){
        ofstream out_u_100;
        unif_metropolis->initialize_acceptance(unif_metropolis->WF100.get(), in_coordinates, wf);  // Initialize Metropolis for 100
        out_u_100.open("./OUTPUT/raggio_unif_w100_pos_"+oss_x.str()+"_"+oss_y.str()+"_"+oss_z.str()+"_"+oss_step.str()+"_step.dat");

        for (int i = 0; i < 100000; i++){ //equilibration
          unif_metropolis->metropolis_method(unif_metropolis->WF100.get(), in_coordinates);
        }
        // Loop over blocks
        for (int j = 0; j < N_block; j++){
          for (int i = 0; i < step_block; i++) {
            unif_metropolis->metropolis_method(unif_metropolis->WF100.get(), in_coordinates);  // Execute Metropolis
            if (l == 3){
              unif_metropolis->print_coordinate_accepted(wf,in_coordinates[0]); //save coordinate
            }
            vec radius = unif_metropolis->get_coord_old();  // Get current coordinates
            r[j] += radius[0];  // Add radius
          }
          r[j] /= step_block;       // Calculate mean radius per block
          r_2[j] = pow(r[j],2);     // Calculate squared mean radius
        }
        // Perform block averaging and write results
        Bloch_mean_method(r, r_2, N_block, out_u_100);
      }

      // If quantum number is 210
      else if( wf == 210){
        ofstream out_u_210;
        unif_metropolis->set_coord_new(in_coordinates);
        unif_metropolis->initialize_acceptance(unif_metropolis->WF210.get(), in_coordinates, wf);  // Initialize Metropolis for 210
        out_u_210.open("./OUTPUT/raggio_unif_w210_pos_"+oss_x.str()+"_"+oss_y.str()+"_"+oss_z.str()+"_"+oss_step.str()+"_step.dat");

        for (int i = 0; i < 100000; i++){
          unif_metropolis->metropolis_method(unif_metropolis->WF100.get(), in_coordinates);
        }
        // Loop over blocks
        for (int j = 0; j < N_block; j++){
          for (int i = 0; i < step_block; i++) {
            unif_metropolis->metropolis_method(unif_metropolis->WF210.get(), in_coordinates);  // Execute Metropolis
            if (l == 3){
              unif_metropolis->print_coordinate_accepted(wf,in_coordinates[0]); //save coordinate
            }
            vec radius = unif_metropolis->get_coord_old();  // Get current coordinates
            r[j] += radius[0];  // Add radius
          }
          r[j] /= step_block;       // Calculate mean radius per block
          r_2[j] = pow(r[j],2);     // Calculate squared mean radius
        }
        // Perform block averaging and write results
        Bloch_mean_method(r, r_2, N_block, out_u_210);
      }
    }
  }
  return 0;
}
