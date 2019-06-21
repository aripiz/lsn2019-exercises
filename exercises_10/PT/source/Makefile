CC = mpicxx
CFLAGS = -Wall -O3 --std=c++14

pt.exe : pt.o random.o tsp.o annealing.o
	$(CC) random.o tsp.o annealing.o pt.o -o ../pt.exe
pt.o : pt.cpp pt.hpp
	$(CC) -c pt.cpp -o pt.o $(CFLAGS)
annealing.o : annealing.cpp annealing.hpp
	$(CC) -c annealing.cpp -o annealing.o $(CFLAGS)
tsp.o : tsp.cpp tsp.hpp
	$(CC) -c tsp.cpp -o tsp.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o ../pt.exe
