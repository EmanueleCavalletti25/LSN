CC = mpic++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o city.o tsp.o
	$(CC) random.o tsp.o city.o main.o -o main.exe 
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
city.o : city.cpp city.h
	$(CC) -c city.cpp -o city.o $(CFLAGS)
tsp.o : tsp.cpp tsp.h
	$(CC) -c tsp.cpp -o tsp.o $(CFLAGS)
clean:
	rm *.o *.dat main.exe seed.out
cleandat:
	rm ./OUTPUT/*.*
esegui: 
	mpiexec -np 5  main.exe