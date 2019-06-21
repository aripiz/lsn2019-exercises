/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019

Genetic algorithm
population.cpp

Francesco Ariele Piziali
*****************************************************************
*****************************************************************/

#include "population.hpp"
#include <algorithm>
#include <cmath>

using namespace std;

Population :: Population(Random& generator, std::vector<City> geometry, int size, int elite_individual, std::vector<double> rates){

  //Random generator
  rnd = &generator;

  //Save realization of the problem
  realization = geometry;
  ngene = realization.size();

  //Save parameters
  npop = size;
  nelite =  elite_individual;
  crossover = rates[0];
  permutation = rates[1];
  tshift = rates[2];
  inversion = rates[3];

  //Create Population
  Populate();
  Sort();

};

void Population :: Populate(){//Create the initial population: a vector of  random individuals
  individual.resize(npop);
  for(int i=0; i<npop; ++i) individual[i] = Route( realization, RandomGenotype() );
}

std::vector<int> Population :: RandomGenotype(){//Create an individual: a route that visits every city once
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

Route Population :: Individual(int i){//Return individual at position i
  return individual[i];
}

double Population :: AveLength(){//Average Length in the population
  double mean = 0;
  for(int i=0; i<npop; ++i){
    mean += individual[i].Length();
  }
  return mean/npop;
}

double Population :: StdLength(){//Standard deviation of the lenght in the population
  double ave = AveLength();
  double sum = 0;
  for(int i=0; i<npop; ++i){
    sum += (individual[i].Length()-ave)*(individual[i].Length()-ave);
  }
  return sqrt(ave/(npop-1));
}

double Population :: MinLength(){//Minimum Length in the population
  return individual.front().Length();
}

double Population :: MaxLength(){//Maximum Length in the population
  return individual.back().Length();
}

int Population :: RandomIndividual(){//Return random individual in the population
  return rnd->Integer(0,npop-1);
}

int Population :: RandomGene(){//Return random gene (city)
  return rnd->Integer(0,ngene-1);
}

void Population :: Breed(){//Reproducion of the population
  vector<Route> children;

  int parent1, parent2;

  for(int i=0; i<npop; ++i){
    parent1 = i;
    do{
      parent2 = RandomIndividual();
    } while(parent2 == i);
    if (rnd->Rannyu()<crossover){
      Crossover(individual[parent1], individual[parent2]);
      children.push_back(individual[parent1]);
    }
    else children.push_back(individual[parent1]);
  }
  individual = children;
}

void Population :: Selection(){//Select individuals that will partecipate in reproduction
  vector<Route> selection;
  double fit_sum = 0;

  for (int i=nelite; i<npop; ++i) fit_sum += individual[i].Fitness();
  for (int i=nelite; i<npop; ++i) individual[i].SetSelectionProbability(fit_sum);

  for (int i=0; i<nelite; ++i) selection.push_back(individual[i]);

  int accepted = nelite;
  while(accepted < npop){
    double p = rnd->Rannyu();
    for(int i=nelite; i<npop; ++i){
      if(accepted >= npop) break;
      if( p < individual[i].SelectionProbability() ) {
        selection.push_back(individual[i]);
        accepted++;
      }
    }
  }

  individual = selection;
}

void Population :: Mutation(){//Mutation of the individuals

  //Permutation of elements at p1 and p2
  int p1 = RandomGene();
  int p2 = RandomGene();
  for(int i=0; i<npop; ++i){
    if(rnd->Rannyu()<permutation) individual[i].Permutation(p1,p2);
  }
  //Total shift of places
  int places = RandomGene();
  for(int i=0; i<npop; ++i){
    if(rnd->Rannyu()<tshift) individual[i].Shift(places);
  }

  //Inversion
  p1 = RandomGene();
  p2 = RandomGene();
  for(int i=0; i<npop; ++i){
    if(rnd->Rannyu()<inversion) individual[i].Inversion(p1,p2);
  }

}

void Population :: Crossover(Route& parent1, Route& parent2){//Crossover of two individuals
  int cut = rnd->Integer(1,ngene-2);
  vector<int> child;
  vector<int> missing;
  for(int i=0; i<cut; ++i) child.push_back(parent1.Genotype()[i]);

  for(int i=cut; i<ngene; ++i){
      for(int j=0; j<ngene; ++j) if(parent2.Genotype()[j] == parent1.Genotype()[i]) missing.push_back(j);
  }
  sort(missing.begin(), missing.end());

  for(int i=0; i<missing.size(); ++i) child.push_back(parent2.Genotype()[missing[i]]);
  parent1.SetGenotype(child);
}

void Population :: Update(){//Update every individual
  for(int i=0; i<npop; ++i) individual[i].Update();
}

void Population :: Print(int igen){
    individual.front().Print(igen);
}

void Population :: Evolve(){//Evolution of the population
  Selection();
  Breed();
  Mutation();
  Update();
  Sort();
}
