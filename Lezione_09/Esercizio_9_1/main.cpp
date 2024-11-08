#include "random.h"
#include "tsp.h"
#include "city.h"

using namespace std;
using namespace arma;

int main(){
    

    tsp TSP;
    
    TSP.initialize();
    TSP.initialize_prop();
    int counter = 0;
    int generation = 350;
    double media_loss;
    double alpha =  15.899;
    for (int i = 0; i< generation; i++){
    
    mat pop = TSP.get_mat_pop();
    mat pop_new = pop;
    arma::vec cum_fit = TSP.cumulative_fit(pop,alpha); //Fitness vector - cumulative; necessary for the select
    
    while (counter < 100){
    int id_1 = 0;
    int id_2 = 0;
    
    id_1 = TSP.select_idx_fit(cum_fit);         //Index vector from cum_fit 
    id_2 = TSP.select_idx_fit(cum_fit);         //Index vector from cum_fit 
        
    vec gen_1 = TSP.get_row_pop(id_1);          // --> generator_1
    vec gen_2= TSP.get_row_pop(id_2);           // --> generator_2
        
    TSP.crossover(gen_1,gen_2);                 //crossover with probability p_c; could be generated crossover son or equal one;
    TSP.mutation(gen_1);                        //mutation (porposed 4) with probability p_m;
    TSP.mutation(gen_2);                        //mutation (porposed 4) with probability p_m;

    pop_new.row(counter) = gen_1.t();
    counter++;
    pop_new.row(counter) = gen_2.t();
    counter++;
    }
    
    counter = 0;
    
    TSP.set_mat_pop(pop_new);
    media_loss=TSP.write_loss_function();
    TSP.write_loss_function_single();
    cout << "generation: " << i+1 << endl;
    TSP.write_best_path(TSP.get_mat_pop());
    
    };
    
    cout << "ideal_potential " << alpha <<endl;
    
    return 0;
}
