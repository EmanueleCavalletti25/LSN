## How to use my code - LSN Exercise 06

### Preparation of data

in `INPUT/ISING/` you have to select:
- in `input_equilibration.dat` and in `input.dat`: you have to select `SIMULATION TYPE ` (2 -->MC, 3 --> GIBBS), you can also exchange other parameters.

### Simulation

The compilation of the simulation itself is slightly different. Go into `Lezione_6/NSL_SIMULATOR/SOURCE` and
- to compile hit `make`
- to remove `*.o` and `*.exe` files, hit `make clean`
- to remove `*.dat` files hit `make remove_ISING`
- finally to execute the simulation, hit `./simulation.exe <method: equilibration/evolution> : ISING` 