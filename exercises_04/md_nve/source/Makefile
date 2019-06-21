CC = g++
CFLAGS = -Wall -O3 --std=c++11

Molecular_Dynamics_NVE.exe : Molecular_Dynamics_NVE.o random.o
	$(CC) random.o Molecular_Dynamics_NVE.o -o Molecular_Dynamics_NVE.exe
Molecular_Dynamics_NVE.o : Molecular_Dynamics_NVE.h Molecular_Dynamics_NVE.cpp
	$(CC) -c Molecular_Dynamics_NVE.cpp -o Molecular_Dynamics_NVE.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Molecular_Dynamics_NVE.exe seed.out
