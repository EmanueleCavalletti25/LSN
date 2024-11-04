//
//  Esercizio_1.cpp
//  
//  Created by Emanuele Cavalletti on 29/02/2024.
//

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "funzioni.hpp"
#include <vector>
#include <cmath>

using namespace std;

int main(int argc, char *argv[]) {
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " <Number_of_throws> <Number_of_blocks>" << endl;
        return -1;
    }

    // Number of throws
    int N = atoi(argv[1]);
    // Number of blocks
    int M = atoi(argv[2]);
    // Number of data per block
    int L = (int)(N / M);

    // Ensure the number of blocks is not greater than the number of throws
    while (M > N) {
        cout << "Number of blocks exceeds the number of generated data." << endl;
        cin >> N;
        cin >> M;
    }

    // Random number generation setup
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } else {
        cerr << "PROBLEM: Unable to open Primes" << endl;
    }
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else {
        cerr << "PROBLEM: Unable to open seed.in" << endl;
    }

    ofstream output_data("Valori_pi.txt");
    if (!output_data.is_open()) {
        cout << "Error opening the file!" << endl;
        return 1; // Exit with error
    }

    // Necessary variables to calculate Pi
    double ago = 1;
    double dist_lines = 1.5;

    // Vectors for calculations
    vector<double> mean_pi;       // Mean
    vector<double> mean_pi_2;     // Mean^2
    vector<double> prog_pi;        // Progressive mean of pi
    vector<double> prog_pi_2;      // Progressive mean of pi^2
    vector<double> error_pi;

    ofstream theta("theta_check.txt");

    double pi = 0;
    double x, y, d;
    int counter;

    // Main loop for blocks
    for (int m = 0; m < M; m++) {
        counter = 0;
        pi = 0;

        for (int i = 0; i < L; i++) {
            // Generate random coordinates within the unit circle
            do {
                x = rnd.Rannyu(); 
                y = rnd.Rannyu(); 
             } while (sqrt(x * x + y * y) > 1);

            theta << atan(x / y) << endl;

            // Randomly generate distance
            d = rnd.Rannyu() * dist_lines;
            if ((d + (ago / 2.) * sin(atan(x / y)) >= dist_lines) || 
                (d - (ago / 2.) * sin(atan(x / y)) <= 0)) {
                counter++;
            }
        }

        // Calculate pi
        pi = (2 * ago * L) / (double)(dist_lines * counter);
        mean_pi.push_back(pi);
        mean_pi_2.push_back(pow(pi, 2));
    }

    // Calculate progressive means
    for (int i = 0; i < M; i++) {
        double prog_sum = 0;
        double prog_sum_2 = 0;

        for (int j = 0; j < i + 1; j++) {
            prog_sum += mean_pi[j];
            prog_sum_2 += mean_pi_2[j];
        }

        // Progressive mean and error calculation
        prog_pi.push_back((prog_sum / (double)(i + 1)));
        prog_pi_2.push_back(prog_sum_2 / (double)(i + 1));
        error_pi.push_back(Error(prog_pi, prog_pi_2, i));
    }

    // Output results
    for (int i = 0; i < mean_pi.size(); i++) {
        output_data << prog_pi[i] << "\t " << error_pi[i] << endl;
        cout << prog_pi[i] << "\t " << error_pi[i] << endl;
    }

    output_data.close();
    theta.close();

    rnd.SaveSeed();
    return 0;
}
