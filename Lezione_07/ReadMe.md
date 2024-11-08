## How to use my code - LSN Exercise 07

### Preparation of data

For the computing of this exercise; I prepare `Input.dat` and `Input_equilibration.dat` file data with the information for each type and method of simulation. To change follow the path example: `INPUT\METHOD\PHASE\input.dat`
where method could be: `MC, MD` (Monte Carlo simulation or Molecular Dynamics simulation) and where phase could be: `GAS,LIQUID,SOLID,ISING`

in `INPUT/MC/ISING/` you have to select:
- in `input_equilibration.dat` and in `input.dat`: you have to select `SIMULATION TYPE ` (2 -->MC, 3 --> GIBBS), you can also exchange other parameters.

### Simulation

The compilation of the simulation itself is slightly different. Go into `Lezione_7/NSL_SIMULATOR/SOURCE` and
- to compile hit `make`
- to remove `*.o` and `*.exe` files, hit `make remove_PHASE_METHOD`
- finally to execute the simulation, hit `./simulation.exe <method: equilibration/evolution> ; PHASE, METHOD` 