## How to use my code - LSN Exercise 05

### Preparation Equilibration

in `NSL_SIMULATOR/INPUT/GAS or LIQUID or SOLID` there are files called
- `input_equilibration.dat` with the information for the equilibration.
- `input.dat` with the information for the equilibration. 
### Simulation


The compilation of the simulation itself is slightly different. Go into `NSL_SIMULATOR/SOURCE` and:
- to compile hit `make`
- to remove `*.o` and `*.exe` files, hit `make clean`
- to remove  files of OUTPUT GAS hit `make remove_GAS`
- to remove  files of OUTPUT LIQUID hit `make remove_LIQUID`
- to remove  files of OUTPUT SOLID hit `make remove_SOLID`
- finally to execute the simulation, hit `./main.exe <method: equilibration or evolution  >`