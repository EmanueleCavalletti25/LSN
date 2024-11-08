#ifndef __Integrale__
#define __Integrale__

#include <cmath>
#include <iostream>
#include "FunzioneBase.h"
#include "random.h"

// Base class for integral calculation using various methods
class Integral {
public:

    // Constructors to initialize integration limits and precision
    Integral(double a, double b) {
        IntegSign(a, b);
        mnstep = 0;
        mh = 0;
        msum = 0;
        m_integral = 0;
        merrore = 0;
        mprec = 0;
    }

    Integral(double a, double b, double prec) {
        IntegSign(a, b);
        mnstep = 0;
        mprec = prec;
        mh = 0;
        msum = 0;
        m_integral = 0;
        merrore = prec + 1;
    }

    // Destructor
    ~Integral() { ; }

    // Returns step size for integration
    double Geth(int N) { return((mb - ma) / (double)N); }
    // Returns the computed integral value
    double GetIntegral() { return m_integral; }

    // Abstract method for integration to be implemented in derived classes
    virtual double Integra(unsigned nstep, FunzioneBase &f, Random &rand) = 0;

protected:

    // Sets integration limits and direction based on input range
    void IntegSign(double a, double b) {
        if (a > b) {
            ma = b;
            mb = a;
            m_sign = -1;
        } else if (a < b) {
            ma = a;
            mb = b;
            m_sign = 1;
        }
        if (a == b) {
            cout << "The integral is zero" << endl;
        }
    }

    // Member variables for integration parameters
    double ma, mb, mh, m_integral, msum, mprec, merrore;
    int m_sign;
    unsigned int mnstep;
};



// Derived class: Uniform sampling method for integration
class Uniforme : public Integral {
public:
    Uniforme(double a, double b) : Integral(a, b) { ; }
    ~Uniforme() { ; }

    virtual double Integra(unsigned int nstep, FunzioneBase &f, Random &rand) override {
        if (nstep < 1) {
            cout << "Negative step size" << endl;
            exit(-1);
        }
        double sum = 0;
        double x;

        for (unsigned int i = 0; i < nstep; i++) {
            x = rand.Rannyu();
            sum += f.Eval(x);  // Adds function value at each random x
        }

        // Calculates the integral using average function value
        m_integral = sum / nstep;
        return m_integral;
    }
};

// Derived class: Importance sampling method for integration
class cos_important_sampling : public Integral {
public:
    cos_important_sampling(double a, double b) : Integral(a, b) { ; }
    ~cos_important_sampling() { ; }

    virtual double Integra(unsigned int nstep, FunzioneBase &f, Random &rand) override {
        if (nstep < 1) {
            cout << "Negative step size" << endl;
            exit(-1);
        }
        double sum = 0;
        double x;

        for (unsigned int i = 0; i < nstep; i++) {
            x = sqrt(1 - rand.Rannyu()) + 1;  // Transforms for importance sampling
            sum += f.Eval(x) / (-2 * x + 2);  // Adjusts by probability density
        }

        // Calculates the integral using importance sampling adjustment
        m_integral = sum / nstep;
        return m_integral;
    }
};

// Derived class: Accept/Reject method for integration
class Acc_Rej : public Integral {
public:
    Acc_Rej(double a, double b) : Integral(a, b) { ; }
    ~Acc_Rej() { ; }

    virtual double Integra(unsigned int nstep, FunzioneBase &f, Random &rand) override {
        if (nstep < 1) {
            cout << "Negative step size" << endl;
            exit(-1);
        }
        double x, r, p;
        double N_hit = 0;
        double p_max = f.Val_max_ass();  // Max probability for accept/reject

        for (unsigned int i = 0; i < nstep; i++) {
            x = rand.Rannyu(ma, mb);
            r = rand.Rannyu();
            p = f.Eval(x);

            if (r < (p / p_max)) {
                N_hit++;
            }
        }

        // Calculates the integral using hit count
        m_integral = m_sign * (N_hit / nstep) * (mb - ma) * p_max;
        return m_integral;
    }
};

#endif // __Integrale__
