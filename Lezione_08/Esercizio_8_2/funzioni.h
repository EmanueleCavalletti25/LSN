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

double Error (const vector<double>& prog_media, const vector<double>& prog_media_2, int n);

//funzioni per fare le medie progressive
vector<double> progressive_method (const vector<double>& media, int blocchi);

void Bloch_mean_method (const vector<double>& media, const vector<double>& media_2, int blocchi, double& mean, double& err);

double distanza_3D (double x, double y, double z);






#endif //__Functions__
