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

//funzioni per fare le medie progressive
vector <double> progressive_method (vector <double> media, int blocchi);

void Bloch_mean_method (vector <double> media, vector <double> media_2 ,int blocchi, ofstream& out );

//costruisco una funzione che mi aiuta a calcolare lo spot price, dato un tempo t
double spot_price_ST_direct (double mean, double volatility, double T, double initial_ST ,Random& rnd);

double spot_price_ST_discreto (double mean, double volatility, double T, double initial_ST ,Random& rnd ,double N_step);


#endif //__Functions__
