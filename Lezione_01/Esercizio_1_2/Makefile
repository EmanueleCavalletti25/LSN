CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o funzioni.o
	$(CC) random.o funzioni.o main.o -o main.exe 
main.o : main.cpp FunzioneBase.h integrali.hpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.hpp
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out

