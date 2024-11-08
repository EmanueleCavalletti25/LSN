//
//  FunzioneBase.h
//  
//
//  Created by Emanuele Cavalletti on 13/09/22.
//

#ifndef FunzioneBase_h
#define FunzioneBase_h
#include <cmath>
#include <iostream>

using namespace std;

const double PI = 3.1415926;
class FunzioneBase {
public:
  
virtual double Eval(double x) const = 0;
virtual double Val_max_ass() const = 0; //ritorna il valore assoluto della funzione in generale

    
virtual ~FunzioneBase(){;};
};

class Xtan : public FunzioneBase {

public:


    Xtan() {a=1;b=1;c=1;};
    Xtan(double x,double y,double z) {a=x; b=y;c=z;};
virtual ~Xtan() {;};
virtual double Eval(double x) const { return x*atan(a*tan(b*(1/x)))-c; }

void SetA (double x)  {a=x;}
void SetB (double y)  {b=y;}
void SetC (double z)  {c=z;}

double GetA() const {return a;}
double GetB() const {return b;}
double GetC() const {return c;}

private:
double a,b,c;
};

class Coseno : public FunzioneBase {

public:


    Coseno() {a=1;b=1;c=1;};
    Coseno(double x,double y,double z) {a=x; b=y;c=z;};
virtual ~Coseno() {;};
virtual double Eval(double x) const { return a*cos(b*x)+c; }
virtual double Val_max_ass() const {return a ;} //ritorna il valore assoluto della funzione in generale

void SetA (double x)  {a=x;}
void SetB (double y)  {b=y;}
void SetC (double z)  {c=z;}

double GetA() const {return a;}
double GetB() const {return b;}
double GetC() const {return c;}

private:
double a,b,c;
};



#endif /* FunzioneBase_h */
