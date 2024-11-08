#ifndef __Functions__
#define __Functions__
#include "random.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <armadillo>


using namespace std;
using namespace arma;

double Error (const vec& prog_media, const vec& prog_media_2, int n);

//funzioni per fare le medie progressive
vec progressive_method (const vec& media, int blocchi);

void Bloch_mean_method (const vec& media, const vec& media_2, int blocchi, ofstream& out );

double distanza_3D (double x, double y, double z);






#endif //__Functions__
