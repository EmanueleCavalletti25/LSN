## How to use my code - LSN Exercise 08

### Preparation of data

For the computing of this exercise; there are two directories with path rispectively `Lezione_8/Esercizio_8_1` and `Lezione_8/Esercizio_8_2`.
In both the `main.cpp` it is possible to set the  `N_block` (number of block) and  `N_step` inside the block.
For `Esercizio_8_1\main.cpp` is possible to set manually the value of `mu` and `sigma`.

### Simulation

The compilation of the simulation itself is slightly different. Go into `Lezione_8/Esercizio_8_*` and
- to compile hit `make`
- to remove `*.o` and `*.exe` files, hit `make clean`
- to remove `OUTPUT` files data, hit `make cleandat`
- finally to execute the simulation, hit `./main.exe <initial_position>`
