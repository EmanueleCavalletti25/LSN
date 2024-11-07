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
      return sqrt((media_2[n] - pow(media[n],2))/n);
   }
};

vector <double> progressive_method(vector <double> media, int blocchi){
   vector <double> prog_media;
    
    for (int i=0; i < blocchi ; i++) {
      
      double prog_sum = 0;
    for (int j = 0; j < i+1; j++){
        
        prog_sum += media[j];
    }
   
   //mean progressive vector
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

   
      out << prog_media[media.size()-1] << "\t " << errore_med[media.size()-1] << endl;
      cout << prog_media[media.size()-1] << "\t " << errore_med[media.size()-1] << " " <<media.size() <<endl;
   media.clear();
   media_2.clear();
   errore_med.clear();
};

double distanza_3D (double x, double y, double z){
   return (x*x + y*y + z*z);
}



   
   
