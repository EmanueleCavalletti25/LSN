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
#include <algorithm>

using namespace std;

int main (int argc, char *argv[]){
   
   // Check if the required arguments are provided
   if (argc < 4){
      cout << "Usage: " << argv[0] << " <Number_Of_Random_Number> <Number_of_blocks> <number of steps>" << endl;
      return -1;
   }
   
   // Number of random numbers to generate
   unsigned int N = atoi(argv[1]);
   // Number of blocks
   unsigned int M = atoi(argv[2]);
   // Data per block
   unsigned int L = (int)(N / M);
   // Number of steps (used for asset price calculation in discrete case)
   unsigned int step = atoi(argv[3]);
   
   // Ensure the number of blocks is not greater than generated data
   while (M > N){
      cout << "Number of blocks exceeds generated data" << endl;
      cin >> N;
      cin >> M;
   }
   
   /////////////////////////////////////////////////////////////////////////////////////////////
   // Random number generator setup
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
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
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
   
   /////////////////////////////////////////////////////////////////
   // Exercise 3.1.1 - Direct Sampling
   /////////////////////////////////////////////////////////////////
   
   // Open output files
   ofstream out_d_C("GBW_direct_Call.txt");
   if (!out_d_C.is_open()) {
         cout << "Error opening file!" << endl;
         return 1; // Exit with error
     }
   ofstream out_d_P("GBW_direct_Put.txt");
   
   // Initial financial data setup
   double S_init = 100.;          // Initial asset price
   double r = 0.1;                // Risk-free interest rate
   double K = 100.;               // Strike price
   double sigma = 0.25;           // Volatility
   double T = 1;                  // Expiration time
   double S = 0;                  // Spot price

   // Vectors to store call and put means and their squares
   vector<double> call_media, call_media_2;
   vector<double> put_media, put_media_2;

   double C = 0;  // Sum for call option
   double P = 0;  // Sum for put option

   // Direct sampling: Loop over blocks and data points
   for (int m = 0; m < M; m++){
      for (int l = 0; l < L; l++){
         S = spot_price_ST_direct(r, sigma, T, S_init, rnd);
         C += exp(-1 * r * T) * max(0., S - K);
         P += exp(-1 * r * T) * max(0., K - S);
      }
      call_media.push_back(C / (double)L);
      call_media_2.push_back(pow(C / (double)L, 2));
      put_media.push_back(P / (double)L);
      put_media_2.push_back(pow(P / (double)L, 2));

      C = 0;
      P = 0;
   }

   // Calculate and write results for direct sample calls and puts
   cout << endl << "GBW DIRECT CALL" << endl;
   Bloch_mean_method(call_media, call_media_2, M, out_d_C); // Block method for Call (direct)
   cout << endl << "GBW DIRECT PUT" << endl;
   Bloch_mean_method(put_media, put_media_2, M, out_d_P); // Block method for Put (direct)

   out_d_C.close();
   out_d_P.close();

   /////////////////////////////////////////////////////////////////
   // Discretized Sampling
   /////////////////////////////////////////////////////////////////

   // Vectors to store means for discretized call and put options
   vector<double> call_media_dis, call_media_2_dis;
   vector<double> put_media_dis, put_media_2_dis;

   ofstream out_dis_C("GBW_discrete_Call.txt");
   ofstream out_dis_P("GBW_discrete_Put.txt");
   
   double C_d = 0;  // Sum for call option (discrete)
   double P_d = 0;  // Sum for put option (discrete)

   // Discretized sampling: Loop over blocks and data points
   for (int m = 0; m < M; m++){
      for (int l = 0; l < L; l++){
         S = spot_price_ST_discreto(r, sigma, T, S_init, rnd, step);
         C_d += exp(-1 * r * T) * max(0., S - K);
         P_d += exp(-1 * r * T) * max(0., K - S);
      }
      call_media_dis.push_back(C_d / (double)L);
      call_media_2_dis.push_back(pow(C_d / (double)L, 2));
      put_media_dis.push_back(P_d / (double)L);
      put_media_2_dis.push_back(pow(P_d / (double)L, 2));

      C_d = 0;
      P_d = 0;
   }

   // Calculate and write results for discrete sample calls and puts
   cout << endl << "GBW DISCRETE CALL" << endl;
   Bloch_mean_method(call_media_dis, call_media_2_dis, M, out_dis_C); // Block method for Call (discrete)
   cout << endl << "GBW DISCRETE PUT" << endl;
   Bloch_mean_method(put_media_dis, put_media_2_dis, M, out_dis_P); // Block method for Put (discrete)

   out_dis_C.close();
   out_dis_P.close();

   // Save random generator seed state
   rnd.SaveSeed();
   return 0;
}
