/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019

Exercise 9.1
Genetic algorithm for TSP
ga.hpp

Francesco Ariele Piziali
*****************************************************************
*****************************************************************/
#include "annealing.hpp"
#include "random.h"

//Random numbers generator
Random rnd("seed.in");

//Realization of the Tsp
int ncities;
std::vector<City> cities_list;

//Population
Annealing * tsp;

//Parameters of the simulation
int nstep;
int ncool;
double beta;
double cool;

//Mutation rates
double permutation;
double inversion;
double crossover;
double pshift;

//Simulation
bool geometry;
int print_best;

//FUNCTIONS
void Initialize();
void Realization();
void Averages(int);
void Finalize();
