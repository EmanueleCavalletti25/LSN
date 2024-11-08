//
//  Esercizio_1.cpp
//  
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


int main (int argc, char *argv[]){
   
   if (argc < 2){
      cout << "Uso programma: " << argv [0] << "<Number_Of_Random_Number>" << endl;
      return -1;
   }
   // numeri casuali da generare
   int N = atoi (argv[1]);
   
   //generazione dei numeri casuali
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();
   
   ifstream input("seed.in");
   string property;
   
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
      
   }
   else cerr << "PROBLEM: Unable to open seed.in" << endl;

   double lambda = 1;
   double mu = 0;
   double sigma = 1;


   vector <double> ripetition = {1,2,10,100};

   ofstream File_lin ;
   ofstream File_lor ;
   ofstream File_exp ;

   File_lin.open("Dati_lineare.txt");
   File_lor.open("Dati_lorentziana.txt");
   File_exp.open("Dati_esponenziali.txt");
   
   // faccio girare con un ciclo for il numero di ripetizioni richieste
   for (int n = 0; n < ripetition.size(); n++){
      
   // tengo conto di quanti dati devo creare
      for (int j = 0; j < N; j++ ){
   //  tengo conto di quante ripetizioni faccio
         double lin=0;
         double lor=0;
         double exp = 0;
         
         for (int s = 0; s<ripetition[n]; s++){
   //scrivo i vari dati ottenuti con le diverse distribuzioni richieste: lineare, esponenziale, lorentziana
            lin += rnd.Rannyu();
            exp += rnd.Exp_distribution(lambda);
            lor += rnd.Cauchy_Lorentz(mu, sigma);
            
            }
            
            lin = lin/ripetition[n];
            File_lin << lin << " ";
            
            exp /= ripetition[n];
            File_exp << exp << " ";
           
            lor = lor/ripetition[n];
            File_lor << lor << " ";
         
         }
         File_exp << endl;
         File_lin << endl;
         File_lor << endl;
      }

   
   File_exp.close();
   File_lin.close();
   File_lor.close();
   
   rnd.SaveSeed();

   return 0;

}
   
   
   
   
   

