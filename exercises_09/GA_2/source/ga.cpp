/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019

Exercise 9.1
Genetic algorithm for TSP
ga.cpp

Francesco Ariele Piziali
*****************************************************************
*****************************************************************/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "ga.hpp"

using namespace std;

/*****************************************************************/
int main(){

  Initialize();

  for(int igen=1; igen<=ngen; ++igen){//Evolve population for ngen generations
    population->Evolve();
    Averages(igen);
  }

  Finalize();

  return 0;
}
/*****************************************************************/

void Realization(){ //Realize an instance of the TSP

  if (geometry == 0){ //Cities on the unit circle
    double theta;
    for (int i = 0; i<ncities; ++i){
      theta = rnd.Rannyu(0,2*M_PI);
      cities_list[i] = City(sin(theta), cos(theta));
    }
  }
  else{ //Cities in a unit square centered in the origin
    for (int i = 0; i<ncities; ++i){
      cities_list[i] = City(rnd.Rannyu(-0.5,0.5), rnd.Rannyu(-0.5,0.5));
    }
  }
}

void Initialize(){//Read input, get a realization of TSP and create a population
  ifstream ReadInput;
  ReadInput.open("input.in");

  ReadInput >> geometry;

  ReadInput >> npop;
  ReadInput >> nelite;
  ReadInput >> ngen;

  ReadInput >> crossover;
  ReadInput >> permutation;
  ReadInput >> tshift;
  ReadInput >> inversion;

  ReadInput >> print_best;
  ReadInput.close();

  cities_list.resize(ncities);

  Realization();

  vector<double> rates({crossover, permutation, tshift, inversion});

  population = new Population(rnd, cities_list, npop, nelite, rates);

  Averages(0);

  cout << "Travelling Salesman Problem - Genetic Algorithm" << endl << endl;

  if(geometry==0) cout <<"Geometry: cities on a unit circle" << endl;
  else cout << "Geometry: cities inside a unit square" << endl;
  cout << "Number of cities to visit = " << ncities << endl;

  cout << endl;
  cout << "Population size = " << npop << endl;
  cout << "Number of elite individuals = " << nelite << endl;
  cout << endl;

  cout << "Number of generations = " << ngen << endl << endl;

  cout << "Rates:" << endl;
  cout << "crossover = " << crossover << endl;
  cout <<  "pair permutation = " << permutation << endl;
  cout <<  "total shift = " <<  tshift << endl;
  cout <<  "inversion = " <<  inversion << endl;
  cout << endl;

  cout << "Initial population:" << endl;
  cout << "Average length = " << population->AveLength() << endl;
  cout << "Minimal length = " << population->MinLength() << endl;
  cout << "Maximal length = " << population->MaxLength() << endl << endl;
}

void Averages(int igen){//Print averages, min and max lengths to file; print configuration to file
  const int wd = 12;

  //Print to file
  ofstream Output;
  Output.open("output/output.len.0",ios::app);
  Output  << setw(wd) << igen
          << setw(wd) << population->AveLength()
          << setw(wd) << population->StdLength()
          << setw(wd) << population->MinLength()
          << setw(wd) << population->MaxLength() << endl;
  Output.close();

  //Show during execution
  if((igen%10 == 0) and (igen!=0 ))
    cout << "Generation " << setw(3) << igen << " - Minimal length = " << setw(5) << population->MinLength() << endl;

  if(igen%print_best == 0) population->Print(igen);
}

void Finalize(){//Show results
  cout << endl;
  cout << "Final population:" << endl;
  cout << "Average length = " << population->AveLength() << endl;
  cout << "Minimal length = " << population->MinLength() << endl;
  cout << "Maximal length = " << population->MaxLength() << endl << endl;

  delete population;
}
