CC =  g++
CFLAGS = -Wall -O3 --std=c++14

sa.exe : sa.o random.o tsp.o annealing.o
	$(CC) random.o tsp.o annealing.o sa.o -o ../sa.exe
sa.o : sa.cpp sa.hpp
	$(CC) -c sa.cpp -o sa.o $(CFLAGS)
annealing.o : annealing.cpp annealing.hpp
	$(CC) -c annealing.cpp -o annealing.o $(CFLAGS)
tsp.o : tsp.cpp tsp.hpp
	$(CC) -c tsp.cpp -o tsp.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o sa.exe seed.out
