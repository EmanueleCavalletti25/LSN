//
//  Esercizio_2_1.cpp
//  
//  Created by Emanuele Cavalletti on 29/02/2024.
//

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "funzioni.hpp"
#include "FunzioneBase.h"
#include "integrali.hpp"
#include <vector>
#include <cmath>

using namespace std;

int main (int argc, char *argv[]){
   
   // Checks if the correct number of arguments is provided
   if (argc < 3){
      cout << "Usage: " << argv[0] << " <Number_Of_Random_Number> <Number_of_blocks>" << endl;
      return -1;
   }
   
   // Number of random numbers to generate
   unsigned int N = atoi(argv[1]);
   // Number of blocks for block averaging
   unsigned int M = atoi(argv[2]);
   // Number of data points per block
   unsigned int L = (int)(N/M);
   
   // Ensures block count doesn't exceed total generated data
   while(M > N){
      cout << "Number of blocks exceeds total generated data" << endl;
      cin >> N;
      cin >> M;
   }
   
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
   // Random number generator setup
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2;
   } else cerr << "ERROR: Unable to open Primes" << endl;
   Primes.close();
   
   ifstream input("seed.in");
   string property;
   
   if (input.is_open()){
      while (!input.eof()){
         input >> property;
         if (property == "RANDOMSEED"){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed, p1, p2);
         }
      }
      input.close();
   } else cerr << "ERROR: Unable to open seed.in" << endl;
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
   // Vectors for storing progressive means and squared means
   vector<double> media_Unif;      // Mean for uniform integration
   vector<double> media_2_Unif;    // Squared mean for uniform integration
   vector<double> media_IS;        // Mean for importance sampling
   vector<double> media_2_IS;      // Squared mean for importance sampling
   
   // Functions and integration methods
   Coseno cos(M_PI/2, M_PI/2, 0);
   Uniforme integrale_cos_Unif(0., 1.);
   cos_important_sampling integrale_cos_IS(0., 1.);

   // Output file setup for uniform and importance sampling integrations
   ofstream output_Unif;
   output_Unif.open("Integrale_uniforme.txt");
   if (!output_Unif.is_open()) {
         cout << "Error opening file!" << endl;
         return 1; // Exit with error
   }
   
   ofstream output_IS;
   output_IS.open("Integrale_Importan_Sampling.txt");
   if (!output_IS.is_open()) {
         cout << "Error opening file!" << endl;
         return 1; // Exit with error
   }
   
   // Calculates integral with uniform and importance sampling methods
   for (int i = 0; i < M; i++){
      double result_Unif = integrale_cos_Unif.Integra(L, cos, rnd);
      media_Unif.push_back(result_Unif);
      media_2_Unif.push_back(pow(result_Unif, 2));

      double result_IS = integrale_cos_IS.Integra(L, cos, rnd);
      media_IS.push_back(result_IS);
      media_2_IS.push_back(pow(result_IS, 2));
   }
   
   // Calculates block averages and errors, and writes to files
   Bloch_mean_method(media_Unif, media_2_Unif, M, output_Unif);
   Bloch_mean_method(media_IS, media_2_IS, M, output_IS);
   
   // Closes output files and saves random generator seed
   output_Unif.close();
   output_IS.close();
   rnd.SaveSeed();
   return 0;
}
