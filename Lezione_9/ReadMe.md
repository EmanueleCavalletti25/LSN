## How to use my code - LSN Exercise 09

### Preparation of data

To set up this exercise, create an `input.dat` file with the necessary configuration details for each simulation type and method. Use the following path to locate the file:

`Lezione_9/Esercizio_9_1/INPUT/input.dat`

This file allows you to specify:

- **SHAPE**: Set to `0` for cities arranged in a circle, or `1` for cities arranged in a square.
- **N_CITY**: The number of cities in the simulation.
- **POPULATION**: The size of the population for the genetic algorithm.

To adjust the **number of generations** and the **loss function power** parameter, α, edit the main code. To set the **mutation** and **crossover** probabilities, make changes in `tsp.h`.


### Simulation

The compilation of the simulation itself is slightly different. Go into `Lezione_9/Esercizio_9_1` and
- to compile hit `make`
- to remove `*.o` and `*.exe` files, hit `make clean`
- finally to execute the simulation, hit `./main.exe`

### hint for simulation:
set this parameter in `Esercizio_9_1/main.cpp`to obtain these result:
- double `alpha = 15.899` perfect circumfernce, `generation = 350`; `population (in input.dat) = 100`
- double `alpha = 12.594` circumference semi-perfect, `generation = 350`; `population (in input.dat) = 100`
- double `alpha = 8`; path without intersection, `generation = 350`; `population (in input.dat) = 400`
