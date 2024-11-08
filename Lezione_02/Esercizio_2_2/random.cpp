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
#include "funzioni.hpp"

using namespace std;

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

vector <double> Random :: RW3D_discreto(unsigned int step ){
    //exctract uniformly a number between 0 and 6, the intervall is divided in 6 part, one for each movement in 3D
    //come back the final distance vector
   double x = 0;
   double y = 0;
   double z = 0;

   vector <double> RW;
   
   for (int i = 0; i < step ; i++){
      double s = Rannyu()*6.;
      
      //condizione su x
      if (s < 1.){ x++; }
      if (s>= 1. && s < 2. ){ x--;}
      //condizione su y
      if (s>= 2. && s < 3. ){ y++;}
      if (s>= 3. && s < 4. ){ y--;}
      //condizione su z
      if (s>= 4. && s < 5. ){ z++;}
      if (s>= 5. && s < 6. ){ z--;}
      
      RW.push_back(distanza_3D(x,y,z));
      
    }
 
 return (RW);
   
};

double Random :: random_angle(void){
   
   double x = 0;
   double y = 0;
   do{

      x = Rannyu(-1,1);
      y = Rannyu(-1,1);
      
   } while ((x*x + y*y) < 1);

   if (y >= 0){
      return acos(x/sqrt(x*x+y*y));
   }

   if (y < 0){
      return (2*M_PI-acos(x/sqrt(x*x+y*y)));
   }
   
   cout << "errore" << endl;

   exit(-11);
};

vector <double> Random :: RW3D_continuum (unsigned int step ){

   double x = 0;
   double y = 0;
   double z = 0;
   vector <double> RW;

   for (int i = 0; i < step; i++){
      
      double u = Rannyu(); // Uniformly distributed random number in [0, 1]
      double v = Rannyu(); // Uniformly distributed random number in [0, 1]
      
      double theta = acos(2 * u - 1); // Polar angle
      double phi = 2 * M_PI * v;      // Azimuthal angle
      
      x += sin(theta)*cos(phi);
      y += sin(theta)*sin(phi);
      z += cos(theta);


      /*ofstream out_angle;
      out_angle.open("angle_gen.dat",ios::app);
      out_angle << sin(theta)*cos(phi) << " " << sin(theta)*sin(phi) << " " << cos(theta) << endl;
      out_angle.close();*/
      
   
      RW.push_back(distanza_3D(x,y,z));
      
   }
   
   return (RW);
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
