CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o funzioni.o metropolis.o wavesfunction.o
	$(CC) random.o metropolis.o wavesfunction.o funzioni.o main.o -o main.exe 
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
metropolis.o : metropolis.cpp metropolis.h
	$(CC) -c metropolis.cpp -o metropolis.o $(CFLAGS)
wavesfunction.o : wavesfunction.cpp wavesfunction.h
	$(CC) -c wavesfunction.cpp -o wavesfunction.o $(CFLAGS)

clean :
	rm *.o *.dat main.exe seed.out
cleandat:
	rm ./OUTPUT/*.*