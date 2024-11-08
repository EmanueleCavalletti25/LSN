## How to use my code - LSN Exercise 05

### Preparation of data

in `main.cpp` you have to select:
- `in_coordinates[0]` ––>  // initial position
- `N_block` =         ––>  // Number of blocks
- `step_block`        ––>  // Step in the block
### Simulation

The compilation of the simulation itself is slightly different. Go into `Esercizio_5_1_new/` and
- to compile hit `make`
- to remove `*.o` and `*.exe` files, hit `make clean`
- to remove `*.dat` files hit `make cleandat`
- finally to execute the simulation, hit `./main.exe <method: gaussian/uniform> wave (n,l,m): <100 or 210>` 