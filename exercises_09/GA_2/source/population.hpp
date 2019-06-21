/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019

Genetic algorithm
population.hpp

Francesco Ariele Piziali
*****************************************************************
*****************************************************************/

#include "tsp.hpp"
#include "random.h"

#ifndef __Population__
#define __Population__

class Population {
private:
  //Random generator
  Random * rnd;

  //Parameters
  int npop;
  int nelite;
  int ngene;

  //Vectors
  std::vector<Route> individual;
  std::vector<City> realization;

  //Rates
  double crossover;
  double tshift;
  double inversion;
  double permutation;

  //Creation
  void Populate();
  std::vector<int> RandomGenotype();

  //Genetic actions
  void Mutation();
  void Selection();
  void Breed();
  int RandomGene();
  int RandomIndividual();
  void Crossover(Route& , Route& );

  //Sorting and update
  void Update();
  void Sort(){sort(individual.begin(), individual.end());} //Sort population by decreasing fitness

public:
  //Creator and destructor
  Population(){};
  Population(Random&, std::vector<City>, int, int, std::vector<double>);
  ~Population(){};

  //Evolution
  void Evolve();

  //Results
  double AveLength();
  double MaxLength();
  double MinLength();
  double StdLength();
  void Print(int);
  Route Individual(int);

};


#endif
