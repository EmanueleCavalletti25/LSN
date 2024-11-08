/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "random.h"
#include "funzioni.h"

using namespace std;
//const double a_0 = 0.0529 ; // Bohr radius constant

Random :: Random(){}
// Default constructor, does not perform any action

Random :: ~Random(){}
// Default destructor, does not perform any action

void Random :: SaveSeed(){
   // This function saves the current state of the random number generator to a file "seed.out"
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << "RANDOMSEED	" << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   // This function generates a random number from a Gaussian distribution with given mean and sigma
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   // This function generates a random number in the range [min, max)
   return min+(max-min)*Rannyu();
}

vec Random :: random_newcoord_unif(vec coord_old, double sigma) {
    double r_new;
    vec newcoord = vec(3, fill::zeros);
    //coord_old.print("coordinate vecchie: ");
    //ofstream out_cor;
    //out_cor.open("./OUTPUT/unif_coordinate_generation_xyz.dat", ios::app);
   

    double r = coord_old[0];      // raggio
    double theta = coord_old[1];  // angolo polare(in radianti)
    double phi = coord_old[2];    // angolo azimutale (in radianti)

    // Conversione da coordinate sferiche a cartesiane
    double x = r * sin(theta) * cos(phi);
    double y = r * sin(theta) * sin(phi);
    double z = r * cos(theta);

    // Generazione casuale di nuove coordinate uniformi all'interno di una sfera unitaria
    double r_temp;
    do {
        r_temp = 0.0; // Inizializza r_temp a 0 per ogni nuovo tentativo
        for (int i = 0; i < 3; i++) {
            newcoord[i] = Rannyu(-1., 1.) * sigma;
            r_temp += pow(newcoord[i], 2);
        }
    } while (sqrt(r_temp) > sigma); // La condizione è ora basata su r_temp
    //out_cor << newcoord[0] <<" "<<newcoord[1] << " " << newcoord[2] << endl;
    x += newcoord[0];
    y += newcoord[1];
    z += newcoord[2];

    r_new = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

    // Conversione inversa da cartesiane a sferiche
    coord_old[0] = r_new;               // r
    coord_old[1] = acos(z / r_new);     // theta (angolo zenitale)
    coord_old[2] = atan2(y, x);         // phi (angolo azimutale)

    //out_cor.close();
    // Aggiornamento delle coordinate
    //coord_old.print("coordinate aggiornate: ");
    return coord_old;
}

vec Random :: random_newcoord_gauss (vec coord_old, double sigma){
   vec newcoord = vec(3, fill::zeros);

   //ofstream out_cor;
   //out_cor.open("./OUTPUT/gauss_coordinate_generation_xyz.dat", ios::app);
   

    double r = coord_old[0];      // raggio
    double theta = coord_old[1];  // angolo azimutale (in radianti)
    double phi = coord_old[2];    // angolo zenitale (in radianti)

    // Conversione da coordinate sferiche a cartesiane
    double x = r * sin(theta) * cos(phi);
    double y = r * sin(theta) * sin(phi);
    double z = r * cos(theta);

    // Generazione casuale di nuove coordinate gaussiane
    for (int i = 0; i < 3; i++) {
        newcoord[i] = Gauss(0.0, sigma );  // Genera un valore gaussiano con media 0 e deviazione standard sigma
    }
   // out_cor << newcoord[0] <<" "<<newcoord[1] << " " << newcoord[2] << endl;
    
    x += newcoord[0];
    y += newcoord[1];
    z += newcoord[2];

    double r_new = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

    // Conversione inversa da cartesiane a sferiche
    coord_old[0] = r_new;               // r
    coord_old[1] = acos(z / r_new);     // theta (angolo zenitale)
    coord_old[2] = atan2(y, x);         // phi (angolo azimutale)

    // Aggiornamento delle coordinate
    //out_cor.close();
    return coord_old;
};



double Random :: Rannyu(void){
  // This function generates a random number in the range [0,1)
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  // This function sets the seed and parameters of the random number generator
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}




/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
