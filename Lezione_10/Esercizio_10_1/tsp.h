#pragma once

#include <iostream>
#include "random.h"
#include <cmath>
#include "city.h"
#include <armadillo>
#include "mpi.h"

using namespace std;
using namespace arma;

class tsp{
    private:
    int _shape;         //choose square shape or circonference
    int _ncity;         // number of city
    double _r;             // radium of circonference or side lenght of square
    int _npop;          // how many repliche works in the algorithm
    double p_m = 0.045;  //probability of mutation
    double p_c = 0.81;  //probability of crossover
    
    public:

    tsp() {};           //constructor
    ~tsp();            //destructor
    //FUNCTION TO INITIALIZE
    
    void initialize(int rank);  //function to initialize the system: shape and other parameters
    void initialize_prop(int rank); //creation of repliche in the system
    Random _rnd;

    //FUNCTION OF GENETIC ALGORITHM
    double Fitness_f_1 (vec x);
    vec Loss_vector(mat pop);
    void pair_permutation(vec& x);              //pair permutation of 2 element [1,2,3,4,5]->[1,4,3,2,5];
    void shift_n(vec& x);                       //shift +n of m contiguos element [1,2,3,4,5]->[1,4,5,2,3];
    void swap_m_city(vec& x);                   //swap of m contiguos element with other m (different) [1,2,3,4,5,6]->[1,5,6,4,2,3];
    void swap_back(vec& v);                     //swap back m region  [1,2,3,4,5,6]->[1,5,4,3,2,6];
    vec cumulative_fit(mat popul , double pot);         //create the vector of cumulative probability of fitness (first part of select algorithm)
    int select_idx_fit(vec cum_fit);            //select the index of cumulative probability of fitness (second part of select algorithm)
    void children_generator(vec& a, vec& b);    //generate vector-son with crossover principia;
    
    void mutation(vec& v);                      // algorithm which take a vector and with a probability p_m of applicate some mutation function
    void crossover (vec& v, vec& w);            // crossover operation
    

    //FUNCTION TO WRITE OUTPUT
    double write_loss_function(int rank);
    void write_best_path(mat pop,int rank);
    void write_core_best(int rank,int min_rank, double global_min);
    double write_loss_function_single(int rank);
    

    //FUNCTION TO ELABORATE DATA
    double distance(vec v, vec w);
    void check_vec(vec v);
    void swap(vec& v ,int id_x ,int id_y);
    double Loss_best_media (vec v);             //media of first 50 better Loss value
    int Loss_best_idx (vec v);                  //index of best value of loss/fitness
    int find_sender_idk(int rank, vec rank_n);

    //FUNCTION FOR PARALLELIZATION:
    void migration (int rank ,vec ranks, vec shuf_rank);

    //VARIABLE UTILIZED IN THE PROGRAMM

    field <city> _cityes;       //field where restored the city and permit to change
    mat m_pop;                  //Matrix where restored all the population created for 
    mat m_dist;                 //Matrix where restored all the distance between citys (will be N* N matrix)
    
    
    //FUNCTION GET  SET
    mat get_mat_pop();
    void set_mat_pop(mat pop_new);
    
    void set_row_pop(vec new_row,int idx_row);  //function to modify the row of matrix population 
    vec get_row_pop(int idx_row);               //function to get the row of matrix population 
    
    double get_pc();

    void set_pc(double p_new);

    double get_pm();

    void set_pm(double p_new);

};