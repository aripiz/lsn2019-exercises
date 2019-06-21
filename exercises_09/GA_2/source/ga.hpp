/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019

Exercise 9.1
Genetic algorithm for TSP
ga.hpp

Francesco Ariele Piziali
*****************************************************************
*****************************************************************/
#include "population.hpp"
#include "random.h"

//Random numbers generator
Random rnd("seed.in");

//Realization of the Tsp
const int ncities = 30;
std::vector<City> cities_list;

//Population
Population * population;

//Parameters of the genetic algorithm
int npop;
int ngen;
int nelite;

//Mutation rates
double permutation;
double tshift;
double crossover;
double inversion;

//Simulation
bool geometry;
int print_best;

//FUNCTIONS
void Initialize();
void Realization();
void Averages(int);
void Finalize();
