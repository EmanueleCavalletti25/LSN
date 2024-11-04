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
}