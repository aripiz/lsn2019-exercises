/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019

Simulated annealing
annealing.hpp

Francesco Ariele Piziali
*****************************************************************
*****************************************************************/

#include "tsp.hpp"
#include "random.h"

#ifndef __Annealing__
#define __Annealing__

class Annealing {
private:
  //Random generator
  Random * rnd;

  //Parameters
  int ngene;

  //Vectors
  Route individual;
  std::vector<City> realization;

  //Rates
  double inversion;
  double pshift;
  double permutation;

  //Creation of individual
  std::vector<int> RandomGenotype();

  //Metropolis
  void Mutation(int);
  int RandomGene();
  double beta;

  //Acceptance
  std::vector<double> accepted, attempted;

public:
  //Creator and destructor
  Annealing(){};
  Annealing(Random&, std::vector<City>, std::vector<double>, double);
  ~Annealing(){};

  //Metropolis
  void Move();
  double Boltzmann(double);

  //Temperature
  void SetBeta(double b){beta = b;}
  double Beta(){return beta;}

  //Results
  double Length(){return individual.Length();}

  void Acceptance();
  void Print(int);

};


#endif
