#include "wavesfunction.h"

const double a_0 = 0.0529;

//0.0529; // Bohr radius

Waves_function::Waves_function() {};

Waves_function::~Waves_function() {};

//Implementation of waves function with quantic number: n = 1; l = 0 ; m = 0
// implementation in unit of bohr, and I used the quadratic versione

Waves_H_100 :: Waves_H_100() {};

double Waves_H_100 :: measure(vec coordinates) {

    return pow(a_0, -1.5)/sqrt(M_PI) * exp(-coordinates[0]);

}
// /a_0


//Implementation of waves function with quantic number: n = 2; l = 1 ; m = 0 ----> r' = r/a_0
Waves_H_210 :: Waves_H_210() {};

double Waves_H_210 :: measure(vec coordinates) {

    return pow(a_0, -3./2.) / 8. * sqrt(2. / M_PI) * coordinates[0] * exp(-coordinates[0] / 2.) * cos(coordinates[1]);

}








