/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019

Exercise 10.1
Simulated annealing for TSP
sa.hpp

Francesco Ariele Piziali
*****************************************************************
*****************************************************************/
#include "annealing.hpp"
#include "random.h"

//Random numbers generator
std::shared_ptr<Random> rnd;

//Parallelization
int numtasks;
int taskid;
double min_length;

//Realization of the Tsp
int ncities;
std::vector<City> cities_list;

//Population
std::shared_ptr<Annealing> tsp;

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
void Info();
void Realization();
void Averages(int);
void Finalize();
void Gather(int);
