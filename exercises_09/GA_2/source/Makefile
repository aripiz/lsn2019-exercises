CC = g++
CFLAGS = -Wall -O3 --std=c++11

ga.exe : ga.o random.o tsp.o population.o
	$(CC) random.o tsp.o population.o ga.o -o ../ga.exe
	rm *.o
ga.o : ga.cpp ga.hpp
	$(CC) -c ga.cpp -o ga.o $(CFLAGS)
population.o : population.cpp population.hpp
	$(CC) -c population.cpp -o population.o $(CFLAGS)
tsp.o : tsp.cpp tsp.hpp
	$(CC) -c tsp.cpp -o tsp.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o ga.exe seed.out
