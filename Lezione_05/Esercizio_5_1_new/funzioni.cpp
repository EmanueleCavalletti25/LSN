#include "random.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "funzioni.h"

using namespace std;
using namespace arma;

/*double Error (vector <double> media, vector<double> media_2, int n){ 
   if (n == 0){
      
      return 0;
   }
   else{
      return sqrt((media_2[n] - pow(media[n],2))/n);
   }
};
*/

   vec progressive_method(const vec& media, int blocchi) {
    vec prog_media(blocchi, fill::zeros); // Inizializza un vettore per le medie progressive

    for (int i = 0; i < blocchi; i++) {
        double prog_sum = 0.0;
        
        // Somma progressiva fino al blocco corrente
        for (int j = 0; j <= i; j++) {
            prog_sum += media[j];
        }
        
        // Calcola la media progressiva per il blocco corrente
        prog_media[i] = prog_sum / (double)(i + 1);
    }

    return prog_media; // Restituisce il vettore con le medie progressive

};

double Error(const vec& prog_media, const vec& prog_media_2, int n) {
    if (n == 0) return 0.0;
    return sqrt((prog_media_2[n] - pow(prog_media[n], 2)) / n);
};

void Bloch_mean_method(const vec& media, const vec& media_2, int blocchi, ofstream& out) {
    // Calcola le medie progressive
    vec prog_media = progressive_method(media, blocchi);
    vec prog_media_2 = progressive_method(media_2, blocchi);

    vec errore_med(blocchi, fill::zeros);

    // Calcola l'errore per ciascun blocco
    for (int i = 0; i < blocchi; i++) {
        errore_med[i] = Error(prog_media, prog_media_2, i);
    // Scrivi l'ultima media progressiva e l'errore corrispondente sul file
    //cout << prog_media[i] << "\t" << errore_med[i] << endl;
     out << prog_media[i] << "\t" << errore_med[i] << endl;
    }
    cout << "valore accettato: "<< prog_media[blocchi-1] << "\t"<<"con errore: " << errore_med[blocchi-1] << endl;
}

double distanza_3D (double x, double y, double z){
   return (x*x + y*y + z*z);
}



   
   
