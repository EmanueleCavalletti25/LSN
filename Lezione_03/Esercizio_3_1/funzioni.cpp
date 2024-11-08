#include "random.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "funzioni.hpp"


using namespace std;

double Error (vector <double> media, vector<double> media_2, int n){ 
   if (n == 0){
      
      return 0;
   }
   else{
      return sqrt((media_2[n] - pow(media[n],2))/(double)n);
   }
};

vector <double> progressive_method(vector <double> media, int blocchi){
   vector <double> prog_media;
    
    for (int i=0; i < blocchi ; i++) {
      //per la media progressiva
      double prog_sum = 0;
    for (int j = 0; j < i+1; j++){
        
        prog_sum += media[j];
    }
   
   //media progressiva_vettore media
      prog_media.push_back((prog_sum/(double)(i+1)));
    }

return prog_media;
};

void Bloch_mean_method (vector <double> media, vector <double>  media_2 ,int blocchi, ofstream& out ){
   
   vector <double> prog_media;
   vector <double> prog_media_2;
   vector <double> errore_med;

   prog_media = progressive_method (media, blocchi);
   prog_media_2 = progressive_method (media_2,blocchi);

   for (int i=0; i < blocchi ; i++) {
      errore_med.push_back(Error (prog_media, prog_media_2, i));
   }

   for (int j = 0; j < media.size(); j++){
      out << prog_media[j] << "\t " << errore_med[j] << endl;
      cout << prog_media[j] << "\t " << errore_med[j] << " " <<media.size() <<endl;
   }
   media.clear();
   media_2.clear();
   errore_med.clear();
   
};

double spot_price_ST_direct (double mean, double volatility, double T, double initial_ST,Random& rnd){
   
   double z = rnd.Gauss(0., 1.);
   return (initial_ST* exp( (mean-0.5*pow(volatility,2)) * T+volatility*z*pow(T,0.5)) );
};

double spot_price_ST_discreto (double mean, double volatility, double T, double initial_ST ,Random& rnd ,double step){
   double z;
   double ST = initial_ST;
   double incremento =  T/ (double) step;
   //double t = 0;

   for (int i = 0; i < step ; i ++){
      z = rnd.Gauss(0., 1.);
      ST = ST * exp ( ( mean-0.5*pow(volatility,2) ) * incremento + volatility * z * pow( incremento ,0.5 ))  ;
      
   }
 
 return ST;

};



   
   
