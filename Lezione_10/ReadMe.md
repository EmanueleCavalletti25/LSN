## How to use my code - LSN Exercise 10

### Preparation of data

To configure this exercise, prepare an `input.dat` file that contains the relevant settings for each type and method of simulation. Use one of the following paths to locate the file:

- `Lezione_10/Esercizio_10_1/INPUT/input.dat`: Here, you can set:
  - **SHAPE**: `0` for cities arranged in a circle, or `1` for cities arranged in a square.
  - **N_CITY**: The number of cities in the simulation.
  - **POPULATION**: The population size for the genetic algorithm.

- `Lezione_10/Esercizio_10_2/INPUT/input.dat`: Here, only **POPULATION** needs to be set.

To adjust the **number of generations** and the **loss function power** parameter, α, edit the main code. Mutation and crossover probabilities can be modified in `tsp.h`.


### Simulation

The compilation of the simulation itself is slightly different. Go into `Lezione_10/Esercizio_10_1` and `Lezione_10/Esercizio_10_1`
- to compile hit `make`
- to execute hit `make esegui`
- to remove `*.o` and `*.exe` files, hit `make clean`
- finally to execute the simulation, hit `./main.exe`

### hint for simulation:
set this parameter in `Esercizio_10_1/main.cpp`to obtain these result:
- double `alpha = 15.899` perfect circumfernce, `generation = 350`; `population (in input.dat) = 100`; `N_migr = 40`
- double `alpha = 8`; path without intersection, `generation = 300`; `population (in input.dat) = 400`; `N_migr = 46`

set this parameter in `Esercizio_10_2/main.cpp`to obtain these result:
- double `alpha = 21` perfect circumfernce, `generation = 900`; `population (in input.dat) = 1300`;
