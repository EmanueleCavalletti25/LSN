#ifndef __Functions__
#define __Functions__
#include "random.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include "FunzioneBase.h"
#include "integrali.hpp"

using namespace std;

double Error (const vector<double>& prog_media, const vector<double>& prog_media_2, int n);

//function to do progressive mean and error
vector <double> progressive_method (const vector<double>& media, int blocchi);

void Bloch_mean_method (vector <double> media, vector <double> media_2 ,int blocchi, ofstream& out );


#endif //__Functions__
