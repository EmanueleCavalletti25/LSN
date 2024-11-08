#include "random.h"
#include "tsp.h"
#include "city.h"
#include "mpi.h"
#include <fstream>

using namespace std;
using namespace arma;

int main(int argc, char* argv[]) {
    
    tsp TSP;
    
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get the total number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get the rank of each process
    double global_min;

    vec rank_n(size, fill::zeros);
    vec rank_n_copy(size, fill::zeros); // Necessary for migration
    for (int i = 0; i < size; i++) {
        rank_n[i] = i; // Initialize rank_n with ranks
    }
    if (rank == 0) {
        rank_n_copy = arma::shuffle(rank_n); // Shuffle ranks only in process 0
    }
    MPI_Bcast(rank_n_copy.memptr(), rank_n_copy.n_elem, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Broadcast shuffled ranks to all processes
    int N_Migr = 90;

    TSP.initialize(rank); // Initialize TSP problem for each process
    TSP.initialize_prop(rank); // Additional TSP initialization for properties

    int counter = 0;
    int generation = 900; // Number of generations for the genetic algorithm
    double media_loss = 0.0;
    double pot = 21; // Parameter for fitness calculation

    if (rank == 0) {
        cout << "potenza " << pot << endl;
    }

    for (int i = 0; i < generation; i++) {
        mat pop = TSP.get_mat_pop(); // Get the current population matrix
        mat pop_new = pop; // Initialize new population matrix
        arma::vec cum_fit = TSP.cumulative_fit(pop, pot); // Fitness vector (cumulative), necessary for selection

        while (counter < TSP.get_popul()) {
            int id_1 = TSP.select_idx_fit(cum_fit); // Select individual based on fitness
            int id_2 = TSP.select_idx_fit(cum_fit); // Select another individual

            vec gen_1 = TSP.get_row_pop(id_1); // Get individual 1
            vec gen_2 = TSP.get_row_pop(id_2); // Get individual 2

            TSP.crossover(gen_1, gen_2); // Apply crossover with probability p_c
            TSP.mutation(gen_1); // Apply mutation with probability p_m
            TSP.mutation(gen_2); // Apply mutation with probability p_m

            pop_new.row(counter) = gen_1.t(); // Insert modified gen_1 into new population
            counter++;
            pop_new.row(counter) = gen_2.t(); // Insert modified gen_2 into new population
            counter++;
        }

        counter = 0;
        TSP.set_mat_pop(pop_new); // Update population matrix with new generation
        media_loss = TSP.write_loss_function(rank); // Calculate and save current loss function value

        // MIGRATION
        if (N_Migr != 0 && (i + 1) % N_Migr == 0) { 
            MPI_Barrier(MPI_COMM_WORLD); // Synchronize before migration

            vec L_V = TSP.Loss_vector(TSP.get_mat_pop()); // Get loss values for current population
            int idx = TSP.Loss_best_idx(L_V); // Get index of best individual

            arma::vec temp_vec = (TSP.get_mat_pop()).row(idx).t(); // Get best individual (row vector)
            arma::vec temp_vec_copy = temp_vec; // Copy of the best individual

            MPI_Send(temp_vec.memptr(), temp_vec.n_elem, MPI_DOUBLE, rank_n_copy[rank], rank, MPI_COMM_WORLD); // Send best individual to target rank
            MPI_Barrier(MPI_COMM_WORLD); // Synchronize after sending

            int sender_rank = TSP.find_sender_idk(rank, rank_n_copy); // Find rank that sends data
            MPI_Recv(temp_vec_copy.memptr(), temp_vec.n_elem, MPI_DOUBLE, sender_rank, sender_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Receive best individual from sender rank
            MPI_Barrier(MPI_COMM_WORLD);

            TSP.set_row_pop(temp_vec_copy, idx); // Update population with received individual
            MPI_Barrier(MPI_COMM_WORLD); // Final synchronization after migration
        }

        MPI_Barrier(MPI_COMM_WORLD);
        TSP.write_best_path(TSP.get_mat_pop(), rank); // Save the best path for each process
    }

    vec L_V_n = TSP.Loss_vector(TSP.get_mat_pop()); // Final loss vector after all generations
    media_loss = TSP.Loss_best_media(L_V_n); // Calculate the best average loss

    int min_rank;
    MPI_Reduce(&media_loss, &global_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD); // Find minimum loss among all processes
    double loss_rank[size]; // Array to store the loss values of each rank
    MPI_Barrier(MPI_COMM_WORLD); // Synchronize all processes

    MPI_Allgather(&media_loss, 1, MPI_DOUBLE, loss_rank, 1, MPI_DOUBLE, MPI_COMM_WORLD); // Collect all media_loss values to all processes
    if (rank == 0) {
        cout << "global min is: " << global_min << endl; // Output the global minimum loss
    }

    if (rank == 0) {
        min_rank = -1; // Initialize to an invalid value
        global_min = loss_rank[0]; // Initialize global_min with the first value

        for (int i = 0; i < size; i++) {
            if (loss_rank[i] <= global_min) {
                global_min = loss_rank[i];
                min_rank = i; // Save the corresponding rank
            }
        }
        TSP.write_core_best(0, min_rank, global_min); // Write the core with the best result
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize(); // Finalize MPI

    return 0;
}
