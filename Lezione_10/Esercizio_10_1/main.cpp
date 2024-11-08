#include "random.h"
#include "tsp.h"
#include "city.h"
#include "mpi.h"
#include <fstream>

using namespace std;
using namespace arma;

int main(int argc, char* argv[]){
    
    tsp TSP;
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get the total number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get the rank of each process

    // Setup rank vector for migration
    vec rank_n(size, fill::zeros);
    vec rank_n_copy(size, fill::zeros); // Copy of rank vector, for shuffling in migration
    for (int i = 0; i < size; i++) rank_n[i] = i; 
    if (rank == 0) rank_n_copy = arma::shuffle(rank_n); // Shuffle only in process 0
    MPI_Bcast(rank_n_copy.memptr(), rank_n_copy.n_elem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    int N_Migr = 40; // Migration frequency in generations
    // 46 for circunference
    // 40 for square

    TSP.initialize(rank); // Initialize TSP problem for each process
    TSP.initialize_prop(rank); // Additional TSP properties initialization
    
    int counter = 0;
    int generation = 300; // Number of generations
    double media_loss = 0.0;
    double pot = 8; // Fitness parameter
    
    if (rank == 0) {
        for (int j = 0; j < size; j++) {
            cout << "rank_n_copy: " << rank_n_copy[j] << endl;
        }
    }

   
    
    for (int i = 0; i < generation; i++) {
        mat pop = TSP.get_mat_pop(); // Current population matrix
        mat pop_new = pop;           // New population matrix
        arma::vec cum_fit = TSP.cumulative_fit(pop, pot); // Fitness vector (cumulative), for selection

        while (counter < 100) {
            int id_1 = TSP.select_idx_fit(cum_fit); // Select an individual index based on fitness
            int id_2 = TSP.select_idx_fit(cum_fit); // Select another individual

            vec gen_1 = TSP.get_row_pop(id_1); // Get first selected individual
            vec gen_2 = TSP.get_row_pop(id_2); // Get second selected individual

            TSP.crossover(gen_1, gen_2); // Apply crossover with probability p_c
            TSP.mutation(gen_1);         // Apply mutation with probability p_m
            TSP.mutation(gen_2);         // Apply mutation with probability p_m

            pop_new.row(counter++) = gen_1.t(); // Insert modified gen_1 into new population
            pop_new.row(counter++) = gen_2.t(); // Insert modified gen_2 into new population
        }

        counter = 0;
        TSP.set_mat_pop(pop_new); // Update population matrix with the new generation
        media_loss = TSP.write_loss_function(rank); // Calculate and save current loss
        //TSP.write_loss_function_single(rank);
        
        // MIGRATION: Every N_Migr generations, trigger migration process
        if ((i + 1) % N_Migr == 0) {
            MPI_Barrier(MPI_COMM_WORLD); // Synchronize before migration
        
            vec L_V = TSP.Loss_vector(TSP.get_mat_pop()); // Get loss values for current population
            int idx = TSP.Loss_best_idx(L_V); // Get index of the best individual
            
            // Send the best individual to the corresponding rank in rank_n_copy
            arma::vec temp_vec = (TSP.get_mat_pop()).row(idx).t(); 
            arma::vec temp_vec_copy = temp_vec;
             
            MPI_Send(temp_vec.memptr(), temp_vec.n_elem, MPI_DOUBLE, rank_n_copy[rank], rank, MPI_COMM_WORLD);
             MPI_Barrier(MPI_COMM_WORLD);
            
            
            int sender_rank = TSP.find_sender_idk(rank, rank_n_copy); // Get rank to receive from
            
            MPI_Recv(temp_vec_copy.memptr(), temp_vec.n_elem, MPI_DOUBLE, sender_rank, sender_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Barrier(MPI_COMM_WORLD);
            TSP.set_row_pop(temp_vec_copy, idx); // Update best row with received row
            MPI_Barrier(MPI_COMM_WORLD); // Synchronize after migration
        }

        TSP.write_best_path(TSP.get_mat_pop(), rank); // Save the best path for each process
        TSP.write_loss_function_single(rank);
    }

    // FINAL LOSS EVALUATION
    vec L_V_n = TSP.Loss_vector(TSP.get_mat_pop()); // Final loss vector after generations
    media_loss = TSP.write_loss_function_single(rank); // Calculate the best average loss
    double global_min;
    int min_rank;

    // Find the global minimum loss among all processes
    MPI_Reduce(&media_loss, &global_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    double loss_rank[size]; // Array to store loss values of each rank

    MPI_Allgather(&media_loss, 1, MPI_DOUBLE, loss_rank, 1, MPI_DOUBLE, MPI_COMM_WORLD); // Collect all media_loss values
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) cout << "global min: " << global_min << endl;
    
    // Determine the rank with the minimum loss
    if (rank == 0) {
        min_rank = -1;
        global_min = loss_rank[0];

        for (int i = 0; i < size; i++) {
            if (loss_rank[i] <= global_min) {
                global_min = loss_rank[i];
                min_rank = i; // Save the corresponding rank with the minimum loss
            }
        }
        TSP.write_core_best(0, min_rank, global_min); // Write the core with the best result
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize(); // Finalize MPI

    return 0;
}


//COMMENT:
//double pot = 15.899 circle_perfect (350)
//double pot = 7.291; city in a square
//double pot = 8; MA CON POPOLAZIONE 400 SQUARE (300 gen)
//double pot = 12.594 , circle semi-perfect