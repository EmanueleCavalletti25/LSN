//
//  Esercizio_1.cpp
//  
//  Created by Emanuele Cavalletti on 29/02/2024.
//

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "funzioni.hpp"
#include <vector>
#include <cmath>

using namespace std;

int main (int argc, char *argv[]) {
   
   // Check for the correct number of command-line arguments
   if (argc < 4) {
      cout << "Usage: " << argv[0] << " <Number_Of_Random_Numbers> <Number_of_Blocks> <Number_of_Steps>" << endl;
      return -1;
   }
   
   // Set up parameters
   unsigned int N = atoi(argv[1]); // Number of random walks to generate
   unsigned int M = atoi(argv[2]); // Number of blocks
   unsigned int L = (int)(N / M);   // Number of data points per block
   unsigned int step = atoi(argv[3]); // Number of steps
   
   // Ensure that number of blocks does not exceed number of generated numbers
   while (M > N) {
      cout << "Number of blocks exceeds the number of generated data" << endl;
      cin >> N;
      cin >> M;
   }
   
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
   // Random number generation
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes"); // Open file containing prime numbers
   if (Primes.is_open()) {
      Primes >> p1 >> p2;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();
   
   // Read random seed from file
   ifstream input("seed.in");
   string property;
   if (input.is_open()) {
      while (!input.eof()) {
         input >> property;
         if (property == "RANDOMSEED") {
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed, p1, p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
   
   ///////////////////////////////////////////////////////////////// Esercizio 2.2.1 /////////////////////////////////////////////////////////////////////////////
   
   // Open output file for discrete random walk data
   ofstream out_d("RW3D_discreto.txt");
   if (!out_d.is_open()) {
      cout << "Error opening file!" << endl;
      return 1; // Exit with error
   }
   
   // Discrete random walk with 100 steps, simulating 10^4 values, divided into 10^2 blocks
   vector<vector<double>> RW; // Matrix to store random walk results
   
   // Generate random walk data
   for (int i = 0; i < N; i++) {
      RW.push_back(rnd.RW3D_discreto(step));
   }
   
   vector<double> somme; // Vector to hold sums
   vector<double> somme_2; // Vector for squared sums
   vector<vector<double>> posizioni_f; // Matrix for final positions
   vector<vector<double>> posizioni_f_2; // Matrix for squared positions
   
   double pos = 0;
   
   // Nested loops to calculate block averages
   for (int l = 0; l < M; l++) {
      for (int j = 0; j < step; j++) {
         for (int k = 0; k < L; k++) {
            pos += RW[k + l * L][j]; // Sum positions
         }
         somme.push_back(sqrt(pos / ((double)L))); // Store average distance
         somme_2.push_back(pow(pos / (double)L, 1)); // Store square average
         pos = 0; // Reset position sum
      }
      posizioni_f.push_back(somme); // Store results for this block
      posizioni_f_2.push_back(somme_2);
      somme.clear(); // Clear for next block
      somme_2.clear();
  }

  vector<double> media; // Vector for final averages
  vector<double> media_2;

  cout << "Discrete random walk" << endl;
   
  // Calculate averages across blocks
  for (int j = 0; j < step; j++) {
      for (int i = 0; i < M; i++) {
         media.push_back(posizioni_f[i][j]);
         media_2.push_back(posizioni_f_2[i][j]);
      }
      Bloch_mean_method(media, media_2, M, out_d); // Calculate block means
      media.clear(); // Clear for next step
      media_2.clear();
  }
  
  out_d.close(); // Close output file
   
   ////////////////////////////////////////////////// Esercizio 2.2.2 /////////////////////////////////////////////////////////////
   // Similar process for continuous random walk
   
   // Open output file for continuous random walk data
   ofstream out_c("RW3D_continuo.txt");
   if (!out_c.is_open()) {
      cout << "Error opening file!" << endl;
      return 1; // Exit with error
   }
   
   // Generate continuous random walk data
   vector<vector<double>> RW_c; // Matrix for continuous random walk
   
   for (int i = 0; i < N; i++) {
      RW_c.push_back(rnd.RW3D_continuum(step));
   }
   
   vector<double> somme_c; // Vector for sums
   vector<double> somme_2_c; // Vector for squared sums
   vector<vector<double>> posizioni_f_c; // Matrix for final positions
   vector<vector<double>> posizioni_f_2_c; // Matrix for squared positions
   
   double pos_c = 0;
   
   // Nested loops for block averages
   for (int l = 0; l < M; l++) {
      for (int j = 0; j < step; j++) {
         for (int k = 0; k < L; k++) {
            pos_c += RW_c[k + l * L][j]; // Sum positions
         }
         somme_c.push_back(sqrt(pos_c / ((double)L))); // Store average distance
         somme_2_c.push_back(pow(pos_c / (double)L, 1)); // Store square average
         pos_c = 0; // Reset position sum
      }
      posizioni_f_c.push_back(somme_c); // Store results for this block
      posizioni_f_2_c.push_back(somme_2_c);
      somme_c.clear(); // Clear for next block
      somme_2_c.clear();
  }
   
  vector<double> media_c; // Vector for final averages
  vector<double> media_2_c;
   
  cout << "Continuous random walk" << endl;
  // Calculate averages across blocks
  for (int j = 0; j < step; j++) {
      for (int i = 0; i < M; i++) {
         media_c.push_back(posizioni_f_c[i][j]);
         media_2_c.push_back(posizioni_f_2_c[i][j]);
      }
      Bloch_mean_method(media_c, media_2_c, M, out_c); // Calculate block means
      media_c.clear(); // Clear for next step
      media_2_c.clear();
  }
  
  out_c.close(); // Close output file

  rnd.SaveSeed(); // Save random seed
  return 0; // Exit program
}

