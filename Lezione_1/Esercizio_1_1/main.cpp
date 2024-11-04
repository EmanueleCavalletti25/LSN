//
//  Esercizio_1.cpp
//  
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
        cout << "Usage: " << argv[0] << " <Number_Of_Random_Numbers> <Number_of_Blocks>" << endl;
        return -1;
    }

    // Random numbers to generate
    int N = atoi(argv[1]);
    // Number of blocks
    int M = atoi(argv[2]);
    // Number of data points per block
    int L = (int)(N / M);
    
    // In case the number of generated numbers is less than the requested groups
    while (M > N) {
        cout << "Number of blocks exceeds the number of generated data" << endl;
        cin >> N;
        cin >> M;
    }

    // Generate random numbers
    Random rnd;
    int seed[4];
    int p1, p2;

    // Open the Primes file
    ifstream Primes("Primes");
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } else {
        cerr << "PROBLEM: Unable to open Primes" << endl;
    }
    Primes.close();

    // Open the seed file
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

    ofstream output_dati("RandomData_dati.txt");

    if (!output_dati.is_open()) {
        cout << "Error opening the file!" << endl;
        return 1; // Exit with error
    }

    output_dati << "Mean\t Error_mean\t Std.Dev\t Err_Dev\n" << endl;

    // Calculation of block averages and progressive values

    // For Part 1 Exercise 1
    vector<double> media;        // mean
    vector<double> media_2;     // mean squared
    vector<double> prog_media;   // progressive mean calculation
    vector<double> prog_media_2; // progressive mean squared calculation
    vector<double> errore_med;   // calculated error of the progressive mean

    // For Part 2 Exercise 1
    vector<double> devstd;       // progressive standard deviation
    vector<double> devstd_2;     // calculation of squared standard deviation
    vector<double> prog_dev;     // progressive standard deviation calculation
    vector<double> prog_dev_2;   // progressive squared standard deviation calculation
    vector<double> errore_dev;   // calculated error of the standard deviation

    // Generating and loading means and standard deviations using the blocking method
    for (int i = 0; i < M; i++) {
        double sum = 0;
        double sum_dev = 0;

        for (int j = 0; j < L; j++) {
            sum += rnd.Rannyu();
            sum_dev += pow(rnd.Rannyu() - 0.5, 2);
        }

        // Loading means and standard deviations for each block
        media.push_back(sum / (double)L);
        devstd.push_back(sum_dev / (double)L);

        // Loading means and standard deviations squared for each block
        media_2.push_back(pow(media[i], 2));
        devstd_2.push_back(pow(devstd[i], 2));
    }

    // Calculation of progressive means
    for (int i = 0; i < M; i++) {
        // For the progressive mean
        double prog_sum = 0;
        double prog_sum_2 = 0;

        // For the progressive standard deviation
        double prog_dev_sum = 0;
        double prog_dev_sum_2 = 0;

        for (int j = 0; j < i + 1; j++) {
            prog_sum += media[j];
            prog_sum_2 += media_2[j];
            // Progressive sums of standard deviation
            prog_dev_sum += devstd[j];
            prog_dev_sum_2 += devstd_2[j];
        }

        // Progressive mean vector
        prog_media.push_back((prog_sum / (double)(i + 1)));
        prog_media_2.push_back(prog_sum_2 / (double)(i + 1));
        errore_med.push_back(Error(prog_media, prog_media_2, i));

        // Progressive standard deviation
        prog_dev.push_back((prog_dev_sum / (double)(i + 1)));
        prog_dev_2.push_back(prog_dev_sum_2 / (double)(i + 1));
        errore_dev.push_back(Error(prog_dev, prog_dev_2, i));
    }

    for (int i = 0; i < media.size(); i++) {
        output_dati << prog_media[i] << "\t " << errore_med[i] << "\t " << prog_dev[i] << "\t " << errore_dev[i] << endl;
        // cout << prog_media[i] << "\t " << errore_med[i]<< "\t " << prog_dev[i] << "\t " << errore_dev[i]<<endl;
    }

    output_dati.close();

    // Exercise 1.1.3, calculation of chi-squared

    double spazi = 100;
    double ripetizioni = 4000;
    N = 10000;
    double atteso = N / spazi;

    output_dati.open("TestChi2.txt");
    if (!output_dati.is_open()) {
        cout << "Error opening the file!" << endl;
        return 1; // Exit with error
    }

    vector<double> bins; //creation of an histograms
    bins.push_back(0);

    for (int i = 0; i < spazi - 1; i++) {
        double cost = 1. / spazi;
        bins.push_back(i * 0.01 + cost);
    }

    L = bins.size();
    for (int t = 0; t < ripetizioni; t++) {
        vector<int> counter(L, 0);

        for (int i = 0; i < N; i++) {
            double x = rnd.Rannyu();

            if (x > 0.99) {
                counter[L - 1] += 1;
            } else {
                for (int m = 0; m < L; m++) {
                    if (bins[m] <= x && x < bins[m + 1]) {
                        counter[m] += 1;
                    }
                }
            }
        }

        // Chi-squared test
        double chi_quadro = 0;
        for (int s = 0; s < spazi; s++) {
            chi_quadro += (pow(counter[s] - atteso, 2)) / atteso;
        }

        output_dati << chi_quadro << endl;
    }
    output_dati.close();

    rnd.SaveSeed();
    return 0;
}




