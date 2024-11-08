CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -larmadillo
all: simulator.exe
simulator.exe : system.o NSL_SIMULATOR.o particle.o random.o
	$(CC) system.o particle.o random.o NSL_SIMULATOR.o -o simulator.exe $(AFLAGS)
NSL_SIMULATOR.o : NSL_SIMULATOR.cpp system.h particle.h
	$(CC) -c NSL_SIMULATOR.cpp -o NSL_SIMULATOR.o $(CFLAGS)
system.o : system.cpp system.h random.h particle.h
	$(CC) -c system.cpp -o system.o $(CFLAGS)
particle.o : particle.cpp particle.h
	$(CC) -c particle.cpp -o particle.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o simulator.exe
remove_GAS_MD:
	rm ../OUTPUT/MD_NVE/GAS/EQUILIBRATION/*.*
	rm ../OUTPUT/MD_NVE/GAS/*.*
	rm ../OUTPUT/MD_NVE/GAS/CONFIG/*.*
remove_LIQUID_MD:
	rm ../OUTPUT/MD_NVE/LIQUID/EQUILIBRATION/*.*
	rm ../OUTPUT/MD_NVE/LIQUID/*.*
	rm ../OUTPUT/MD_NVE/LIQUID/CONFIG/*.*
remove_SOLID_MD:
	rm ../OUTPUT/MD_NVE/SOLID/EQUILIBRATION/*.*
	rm ../OUTPUT/MD_NVE/SOLID/*.*
	rm ../OUTPUT/MD_NVE/SOLID/CONFIG/*.*
remove_GAS_MC:
	rm ../OUTPUT/MC_NVT/GAS/EQUILIBRATION/*.*
	rm ../OUTPUT/MC_NVT/GAS/*.*
	rm ../OUTPUT/MC_NVT/GAS/CONFIG/*.*
remove_LIQUID_MC:
	rm ../OUTPUT/MC_NVT/LIQUID/EQUILIBRATION/*.*
	rm ../OUTPUT/MC_NVT/LIQUID/*.*
	rm ../OUTPUT/MC_NVT/LIQUID/CONFIG/*.*
remove_SOLID_MC:
	rm ../OUTPUT/MC_NVT/SOLID/EQUILIBRATION/*.*
	rm ../OUTPUT/MC_NVT/SOLID/*.*
	rm ../OUTPUT/MC_NVT/SOLID/CONFIG/*.*
remove_ISING_MC:
	rm ../OUTPUT/MC_NVT/ISING/*.*
	rm ../OUTPUT/MC_NVT/ISING/CONFIG/*.*
	rm ../OUTPUT/MC_NVT/ISING/EQUILIBRATION/*.*
	
