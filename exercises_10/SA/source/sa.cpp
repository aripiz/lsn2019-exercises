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

#include "sa.hpp"

using namespace std;

/*****************************************************************/
int main(){

  Initialize();

  for(int istep=1; istep <= nstep; ++istep){//Simulation

    if (istep%ncool == 0){//Cooling
      beta = beta / cool;    //pow(0.5, istep/ncool);
      tsp->SetBeta(beta);
    }

    tsp->Move();            //Metropolis move

    Averages(istep);        //Print results
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

  ReadInput >> ncities;
  ReadInput >> geometry;

  ReadInput >> nstep;
  ReadInput >> beta;
  ReadInput >> cool;
  ReadInput >> ncool;

  ReadInput >> permutation;
  ReadInput >> inversion;
  ReadInput >> pshift;

  ReadInput >> print_best;
  ReadInput.close();

  cities_list.resize(ncities);

  Realization();

  vector<double> rates({permutation, inversion, pshift});

  tsp = new Annealing(rnd, cities_list, rates, beta);

  Averages(0);

  cout << "Travelling Salesman Problem - Simulated Annealing" << endl << endl;

  if(geometry==0) cout <<"Geometry: cities on a unit circle" << endl;
  else cout << "Geometry: cities inside a unit square" << endl;
  cout << "Number of cities to visit = " << ncities << endl;

  cout << "Number of Metropolis steps  = " << nstep << endl;
  cout << "Initial beta = " << beta << endl;
  cout << "Cooling rate = " << cool << " every " << ncool << " steps " << endl << endl;

  cout << "Rates:" << endl;
  cout << "pair permutation = " << permutation << endl;
  cout << "inversion = " <<  inversion << endl;
  cout << "partial shift = " <<  pshift << endl;
  cout << endl;

  cout << "Start - Length = " << tsp->Length() << endl << endl;
}

void Averages(int istep){//Print averages, min and max lengths to file; print configuration to file
  const int wd = 12;

  //Print to file
  ofstream Output;
  Output.open("output/output.len.0",ios::app);
  Output  << setw(wd) << istep
          << setw(wd) << tsp->Beta()
          << setw(wd) << tsp->Length()
          << endl;
  Output.close();

  //Show during execution
  int iprint = nstep/10;
  if((istep%iprint == 0) and (istep!=0 )){
    cout  << "Step " << setw(3) << istep
          << " - Beta = " << setw(3) << tsp->Beta()
          << " - Length = " << setw(5) << tsp->Length()
          << endl;

    tsp->Acceptance();
  }

  if(istep%print_best == 0) tsp->Print(istep);
}

void Finalize(){//Show results
  cout << "Final - Length = " << tsp->Length() << endl << endl;

  delete tsp;
}
