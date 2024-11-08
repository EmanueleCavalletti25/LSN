#ifndef __Functions__
#define __Functions__
#include "random.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>


using namespace std;

double Error (vector <double> media ,vector <double>  ,int);

//function for mean progressive vector
vector <double> progressive_method (vector <double> media, int blocchi);

void Bloch_mean_method (vector <double> media, vector <double> media_2 ,int blocchi, ofstream& out );

double distanza_3D (double x, double y, double z);

double random_angle();


#endif //__Functions__
