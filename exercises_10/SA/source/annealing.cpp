/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019

Simulated Annealing
annealing.cpp

Francesco Ariele Piziali
*****************************************************************
*****************************************************************/

#include "annealing.hpp"
#include <algorithm>
#include <iostream>
#include <cmath>

using namespace std;

Annealing :: Annealing(Random& generator, std::vector<City> geometry, std::vector<double> rates, double b ){

  //Random generator
  rnd = &generator;

  //Save realization of the problem
  realization = geometry;
  ngene = realization.size();

  //Save parameters
  permutation = rates[0];
  inversion = rates[1];
  pshift = rates[2];

  //Create individual
  individual = Route(realization, RandomGenotype());
  beta = b;

  //Acceptance
  int nacpt = rates.size() + 1;
  accepted.resize(nacpt);
  attempted.resize(nacpt);
  for(int i=0; i<nacpt; ++i){
    accepted[i] = 0;
    attempted[i] = 0;
  }
  
};

std::vector<int> Annealing :: RandomGenotype(){//Create an individual: a route that visits every city once
  vector<int> chromosome(ngene);

  for(int i=0; i<ngene; ++i) chromosome[i] = i; //increasing order of visiting from 0 to ncities-1

  int p1, p2, temp;

  for(int i=0; i<ngene; ++i){ //random couple permutations on the ordered list
    p1 = rnd->Integer(0,ngene-1);
    p2 = rnd->Integer(0,ngene-1);
    temp = chromosome[p1];
    chromosome[p1] = chromosome[p2];
    chromosome[p2] = temp;
  }

  return chromosome;
}

int Annealing :: RandomGene(){//Return random gene (city)
  return rnd->Integer(0,ngene-1);
}

void Annealing :: Mutation(int type){//Mutation of the individuals

    int p1 = RandomGene();
    int p2 = RandomGene();
    int places = RandomGene();
    if (type == 1){ //Permutation of elements at p1 and p2
      if(rnd->Rannyu()<permutation) individual.Permutation(p1,p2);
    }

    if (type == 2){//Inversion shift of places
      if(rnd->Rannyu()<inversion) individual.Inversion(p1,p2);
    }

    if (type == 3){//Partial shift of places
      if(rnd->Rannyu()<pshift) individual.PartialShift(p1,p2,places);
    }
}

void Annealing :: Move(){

  double p, p_old, p_new;
  Route old = individual;

  int type = rnd->Integer(1,2);

  //Old probability
  p_old = Boltzmann(old.Length());

  //Mutation
  Mutation(type);
  individual.Update();

  //New probability
  p_new = Boltzmann(individual.Length());

  //Metropolis test
  p = p_new/p_old;

  if(rnd->Rannyu() <= p){//accept
    accepted[0] += 1.0;     //global
    accepted[type] += 1.0;  //by mutation type

  }
  else{//reject
    individual = old;
  }

  attempted[0] += 1.0;      //global
  attempted[type] += 1.0;   //by mutation type
}

double Annealing :: Boltzmann(double length){
  return exp(-beta*length);
}

void Annealing :: Print(int istep){
    individual.Print(istep);
}

void Annealing :: Acceptance(){
  cout << "Acceptance rates:" << endl;
  cout << "total = " << accepted[0]/attempted[0] << endl;
  cout << "permutation = " << accepted[1]/attempted[1] << endl;
  cout << "inversion = " << accepted[2]/attempted[2] << endl;
  cout << "partial shift = " << accepted[3]/attempted[3] << endl;
  cout << endl;
}
