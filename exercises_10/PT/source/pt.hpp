/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019

Exercise 10.2
Parallel tempering for TSP
pt.hpp

Francesco Ariele Piziali
*****************************************************************
*****************************************************************/
#include "annealing.hpp"
#include "random.h"

//Random numbers generator
std::shared_ptr<Random> rnd;

//Parallelization
int numtasks, taskid;
double min_length;
double accepted_swap;
double attempted_swap;


//Realization of the Tsp
int ncities;
std::vector<City> cities_list;

//Population
std::shared_ptr<Annealing> tsp;

//Parameters of the simulation
int nstep;
int nbeta;
std::vector<double> beta;
int nswap;

//Mutation rates
double permutation;
double inversion;
double crossover;
double pshift;

//Simulation
bool geometry;

//FUNCTIONS
void Initialize();
void Realization();
void Averages(int);
void Finalize();
void Swap(int);
void Info();
