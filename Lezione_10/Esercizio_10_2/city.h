#pragma once
#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <vector>
#include <sstream>
#include "random.h"

using namespace std;
using namespace arma;

class city {
    private:

    const int _ndim = 2;

    vec _pos; //vector position of city;
    vec _pos_old; // old position of the ciy
    string name;
    
    public:
    city(){};         //constructor
    ~city(){}; //destructor

    void initialize();               //initialize property of the city
    

    vec get_pos(); //return the position of city;
    void set_pos (vec pos_new); //set the position of the city;

    vec get_pos_old(); //return the  old position of city;
    void set_pos_old (vec pos_old); //set the  new old position of the city;

    string get_name(); //return the  name of city;
    void set_name (string name_new); //set the name of the city;
};



