#include "random.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "funzioni.h"

using namespace std;
using namespace arma;

vector<double> progressive_method(const vector<double>& media, int blocchi) {
    vector<double> prog_media(blocchi,0.0); // Inizializza un vettore per le medie progressive

    for (int i = 0; i < blocchi; i++) {
        double prog_sum = 0.0;
        
        // Somma progressiva fino al blocco corrente
        for (int j = 0; j <= i; j++) {
            prog_sum += media[j];
        }
        
        // Calcola la media progressiva per il blocco corrente
        prog_media[i] = prog_sum / (double)(i + 1);
        //cout << "prog_media: " << prog_media[i] << endl;
    }

    return prog_media; // Restituisce il vettore con le medie progressive

};


double Error(const vector<double>& prog_media, const vector<double>& prog_media_2, int n) {
    if (n == 0) return 0.0;
    return sqrt((prog_media_2[n] - pow(prog_media[n], 2)) / n);
};

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

double distanza_3D (double x, double y, double z){
   return (x*x + y*y + z*z);
}



   
   
