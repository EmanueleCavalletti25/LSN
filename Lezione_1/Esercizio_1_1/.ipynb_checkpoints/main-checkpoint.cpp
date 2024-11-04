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
   
   if (argc < 3){
      cout << "Uso programma: " << argv [0] << "<Number_Of_Random_Number><Number_of_blocks>" << endl;
      return -1;
   }
   // numeri casuali da generare
   int N = atoi (argv[1]);
   // numeri di blocchi
   int M = atoi (argv[2]);
   // numero di dati per blocco
   int L = (int)(N/M);
   // nel caso in cui i numeri generati siamo siano minori dei gruppi richiesti
   while(M>N){
      cout <<"numeri di blocchi superiore al numero di dati generati" << endl;
      cin >> N;
      cin >> M;
   }
   
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
   
   
   ofstream output_dati ("RandomData_dati.txt");
   
   if (!output_dati.is_open()) {
         cout << "Errore nell'apertura del file!" << endl;
         return 1; // Uscita con errore
     }
   
   output_dati << "Media\t Errore_med\t Dev.Std\t Err_Dev\n" << endl;
   
   //calcolo delle medie dei blocchi e delle cose progressive
   
   
   //Per Parte 1 esercizio 1
   vector <double> media;        //media
   vector <double> media_2;      //media al quadrato
   vector <double> prog_media;   //calcolo della media progressiva
   vector <double> prog_media_2; //calcolo della medie quadrate progressive
   vector <double> errore_med;   //errore calcolato della media progressiva
   
   //Per Parte 2 esercizio 1
   vector <double> devstd;       //Dev.Standard progressiva
   vector <double> devstd_2;     //calcolo della dev. standard quadrate
   vector <double> prog_dev;     //calcolo della dev.standard progressive
   vector <double> prog_dev_2;   //calcolo della dev.standard quadrate progressive
   vector <double> errore_dev;   //errore calcolato della deviazione standard

    
   //generazione e caricamento delle medie e dev.standard, usando il metodo a blocchi
   for (int i=0; i < M ; i++) {
      double sum = 0;
      double sum_dev = 0;
            for(int j=0; j < L; j++){
               sum += rnd.Rannyu();
               sum_dev += pow(rnd.Rannyu()-0.5,2);
      }
   
   //caricamento medie e devst per ogni blocco
      media.push_back (sum/(double)L);
      devstd.push_back(sum_dev/(double)L);
      
   //caricamento medie e devst quadrate per ogni blocco
      media_2.push_back (pow(media[i],2));
      devstd_2.push_back(pow(devstd[i],2));
   }
   
   //Calcolo delle medie progressive
   for (int i=0; i < M ; i++) {
      //per la media progressiva
      double prog_sum = 0;
      double prog_sum_2 = 0;
      //per la deviazione standard progressiva
      double prog_dev_sum = 0;
      double prog_dev_sum_2= 0;

      for (int j = 0; j < i+1; j++){
          prog_sum += media[j];
          prog_sum_2 += media_2[j];
   //somme progressive della deviazione standard
          prog_dev_sum += devstd[j] ;
          prog_dev_sum_2+= devstd_2[j];
      }
   
   //media progressiva_vettore media
      prog_media.push_back((prog_sum/(double)(i+1)));
      prog_media_2.push_back(prog_sum_2/(double)(i+1));
      errore_med.push_back(Error (prog_media, prog_media_2, i));
   
   //dev.std progressiva
      prog_dev.push_back((prog_dev_sum/(double)(i+1)));
      prog_dev_2.push_back(prog_dev_sum_2/(double)(i+1));
      errore_dev.push_back(Error (prog_dev,prog_dev_2, i));
   
   }
   
   for (int i=0; i< media.size(); i++){
      output_dati << prog_media[i] << "\t " << errore_med[i] << "\t " <<prog_dev[i] << "\t " << errore_dev[i]<<endl;
      //cout << prog_media[i] << "\t " << errore_med[i]<< "\t " <<prog_dev[i] << "\t " << errore_dev[i]<<endl;
   }
   
   output_dati.close();


   //esercizio 1.1.3, calcolo del chi quadro

   double spazi = 100;
   double ripetizioni = 100;
   N = 10000;
   double atteso = N/spazi;


output_dati.open ("TestChi2.txt");
   if (!output_dati.is_open()) {
      cout << "Errore nell'apertura del file!" << endl;
      return 1; // Uscita con errore
     }

vector <double> bins;

bins.push_back(0);

for (int i = 0; i < spazi-1; i++){
   double cost = 1./spazi;
   bins.push_back (i*0.01 + cost);
}



L = bins.size();
   for (int t=0; t < ripetizioni; t++){
      vector <int> counter (L,0);

      for (int i = 0; i < N; i++){
         double x = rnd.Rannyu();

         if (x > 0.99){
            counter[L-1] += 1;
         }
         else {
          for (int m = 0; m < L; m++){
               if (bins[m]<= x && x < bins[m+1]){
               counter[m] += 1;
               }
            }
         }
       }
   
   
   //test del chi quadro
   double chi_quadro = 0;
   for (int s = 0; s < spazi; s++){
      chi_quadro = (pow(counter[s]-atteso,2))/atteso;
   }

   output_dati << chi_quadro << endl;

}
   output_dati.close();

   rnd.SaveSeed();
   return 0;
}




/*
vector <double> Metodo_a_Blocchi (vector <double> a, double n_blocchi){
      
   vector <double> prog_a;
   for (int i=0; i < n_blocchi ; i++) {
      //per la quantità progressiva
      double prog_sum = 0;
      //carico le diverse medie progressive
      for (int j = 0; j < i+1; j++){
         prog_sum += a[j];
      }
      prog_a.push_back((prog_sum/(double)(i+1)));
   }
    return prog_a; 
      
}
 */


/*to do list
raccogliere i dati a blocchi ed analizzarli a blocchi, salvare il risultato di medi a varianza in un array, cioè in un vettore
successivamente salvare tutto da una parte.
rifare il conto progressivo e scrivere quante volte fare questo conto
*/

