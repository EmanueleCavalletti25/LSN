#include "random.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "funzioni.hpp"
#include "FunzioneBase.h"
#include "integrali.hpp"

using namespace std;

double Error(const vector<double>& prog_media, const vector<double>& prog_media_2, int n) {
    if (n == 0) return 0.0;
    return sqrt((prog_media_2[n] - pow(prog_media[n], 2)) /(double) n);
};

vector <double> progressive_method(const vector<double>& media, int blocchi){
   vector<double> prog_media(blocchi,0.0);
    for (int i=0; i < blocchi ; i++) {
      //progressive mean
      double prog_sum = 0.;
    for (int j = 0; j <= i; j++){
        prog_sum += media[j];
    }
   
   //progressive mean
      prog_media[i] = prog_sum / (double)(i + 1);
    }

return prog_media;
};

void Bloch_mean_method (vector <double> media, vector <double> media_2 ,int blocchi, ofstream& out ){
   
   vector <double> prog_media;
   vector <double> prog_media_2;
   vector <double> errore_med;

   prog_media = progressive_method (media, blocchi);
   prog_media_2 = progressive_method (media_2,blocchi);

   for (int i=0; i < blocchi ; i++) {
      errore_med.push_back(Error (prog_media, prog_media_2, i));
   }

   for (int i=0; i< media.size(); i++){
      out << prog_media[i] << "\t " << errore_med[i] << endl;
      cout << prog_media[i] << "\t " << errore_med[i] <<endl;
   }

};

///////////////////////////////////////////////////////////////////////////////////////////////////////




/*

void Bloch_mean_method(const vector<double>& media, const vector<double>& media_2, int blocchi, ofstream& out ) {
    
    // Calcola le medie progressive
    vector<double> prog_media = progressive_method(media , blocchi);
    vector<double> prog_media_2 = progressive_method(media_2, blocchi);

    vector<double> errore_med(blocchi,0.0);

    out << "AVERAGE " << setw(12) << "ERRORE: " << endl;
    // Calcola l'errore per ciascun blocco
    for (int i = 0; i < blocchi; i++) {
        errore_med[i] = Error(prog_media, prog_media_2, i);
    
    out <<prog_media[i] << setw(12) << errore_med[i] << endl;
     
    }

    cout << "valore accettato: "<< prog_media[blocchi-1] << "\t"<<"con errore: " << errore_med[blocchi-1] << endl;
}
*/