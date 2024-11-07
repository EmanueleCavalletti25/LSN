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
   
   if (argc < 4){
      cout << "Uso programma: " << argv [0] << "<Number_Of_Random_Number> <Number_of_blocks> <number of steps>" << endl;
      return -1;
   }
   
   // numeri casuali da generare
   unsigned int N = atoi (argv[1]);
   // numeri di blocchi
   unsigned int M = atoi (argv[2]);
   // numero di dati per blocco
   unsigned int L = (int)(N/M);
   //numero di steps
   unsigned int step = atoi(argv[3]);
   // nel caso in cui i numeri generati siamo siano minori dei gruppi richiesti
   while(M>N){
      cout <<"numeri di blocchi superiore al numero di dati generati" << endl;
      cin >> N;
      cin >> M;
   }
   
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
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
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   /*
   //Esercizio 2.2.1
   
   //apertura file
   ofstream out_d ("RW3D_discreto.txt");
   if (!out_d.is_open()) {
         cout << "Errore nell'apertura del file!" << endl;
         return 1; // Uscita con errore
     }
   
   //Il random walk in tempo discreto si costituituisce di 100 step, per ogni step simulo 10^4 valori, e ogni simulazione si divide in 10^2 blocchi

   vector < vector <double> > RW;
   //genero una matrice che mi permette di sommare i vettori
   
      for (int i = 0; i < N; i++){
       
       RW.push_back (rnd.RW3D_discreto(step));
   }
   
   vector <double> somme;
   vector <double> somme_2;
   vector <vector <double> > posizioni_f;
   vector <vector <double> > posizioni_f_2;
   
   double pos = 0;
   
   // 3 cicli: il primo sui blocchi, il secondo sulle colonne, il terzo sulle righe;
   // riduzione a blocchi
  
  for (int l = 0; l < M ; l ++){
      for (int j = 0; j < step ; j++){
         for (int k = 0; k < L; k++){
            pos += RW[k + l*L][j] ;
         }
         somme.push_back(sqrt(pos/((double)L)));
         somme_2.push_back (pow(pos/(double)L,1));
         pos = 0;
      }
      
      posizioni_f.push_back (somme);
      posizioni_f_2.push_back (somme_2);
      somme.clear();
      somme_2.clear();
  }

  vector <double> media;
  vector <double> media_2;
   
   for (int j= 0; j < step ; j++){
      for (int i = 0; i < M; i++){
         media.push_back(posizioni_f[i][j]);
         media_2.push_back(posizioni_f_2[i][j]);
      }
      
      Bloch_mean_method(media,media_2,M,out_d);
      media.clear();
      media_2.clear();
  }
  
  //lavoro sulle medie a blocchi
      
out_d.close();
   
   
   */


   ////////////////////////////////////////////////// Esercizio 2.2.2 /////////////////////////////////////////////////////////////
   // lavoro molto simile, conti pure, ma siamo nel continuo, per cui di fatti, l'unica cosa da modificare è come viene modificato il passo
//Esercizio 2.2.2
   
   //apertura file
   ofstream out_c ("RW3D_continuo.txt");
   if (!out_c.is_open()) {
         cout << "Errore nell'apertura del file!" << endl;
         return 1; // Uscita con errore
     }
   
   //Il random walk in tempo discreto si costituituisce di 100 step, per ogni step simulo 10^4 valori, e ogni simulazione si divide in 10^2 blocchi

   vector < vector <double> > RW_c;
   //genero una matrice che mi permette di sommare i vettori
   
      for (int i = 0; i < N; i++){
        cout << "sei arrivato quiiiiiiiii" << endl;
       RW_c.push_back (rnd.RW3D_continuum(step));
   }
   
   vector <double> somme_c;
   vector <double> somme_2_c;
   vector <vector <double> > posizioni_f_c;
   vector <vector <double> > posizioni_f_2_c;
   
   double pos_c = 0;
   
   // 3 cicli: il primo sui blocchi, il secondo sulle colonne, il terzo sulle righe;
   // riduzione a blocchi
  
  for (int l = 0; l < M ; l ++){
      for (int j = 0; j < step ; j++){
         for (int k = 0; k < L; k++){
            pos_c += RW_c[k + l*L][j] ;
         }
         somme_c.push_back(sqrt(pos_c/((double)L)));
         somme_2_c.push_back (pow(pos_c/(double)L,1));
         pos_c = 0;
      }
     
      posizioni_f_c.push_back (somme_c);
      posizioni_f_2_c.push_back (somme_2_c);
      somme_c.clear();
      somme_2_c.clear();
  }
   
  vector <double> media_c;
  vector <double> media_2_c;
   
   for (int j= 0; j < step ; j++){
      for (int i = 0; i < M; i++){
         media_c.push_back(posizioni_f_c[i][j]);
         media_2_c.push_back(posizioni_f_2_c[i][j]);
      }
      
      Bloch_mean_method(media_c,media_2_c,M,out_c);
      media_c.clear();
      media_2_c.clear();
  }
  
  //lavoro sulle medie a blocchi
      
   out_c.close();

   rnd.SaveSeed();
   return 0;
   
}


