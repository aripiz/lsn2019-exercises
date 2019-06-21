/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019

Travelling Salesman Problem
tsp.cpp

Francesco Ariele Piziali
*****************************************************************
*****************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "tsp.hpp"

using namespace std;
/*****************************************************************/
//Methods of City, Route

double City :: Distance(City b){//Distance from self to city b
 double dx, dy;
 dx = x - b.X();
 dy = y - b.Y();
 //return sqrt(dx*dx + dy*dy);            //L(1)
 return (dx*dx + dy*dy);    //L(2)
}

void City :: Show(){//Show coordinates of city
  const int wd=12;
  cout << setw(wd) << X() << setw(wd) << Y() << endl;
}

Route :: Route(std::vector<City> realization, std::vector<int> visiting_order){//Creator of Route
 cities = realization;
 genotype = visiting_order;
 ncities = cities.size();

 length = SetLength();
 fitness = Fitness();
}

double Route :: SetLength(){//Set Length of the route
 double path = 0;
 int from, to;
 for (int i = 0; i<ncities; ++i){
   from = genotype[i];
   to = genotype[Pbc(i+1)];
   path += cities[from].Distance(cities[to]);
 }

 return path;
}

int Route :: Pbc(int i) {//Algorithm for periodic boundary conditions
 if(i >= ncities) i = i - ncities;
 //else if(i < 0) i = i + ncities;
 return i;
}

void Route :: Show(){//Show order of visiting
 cout << "[";
 for(int i=0; i<ncities-1; ++i) cout << genotype[i] << " " ;
 cout << genotype.back() << "]" << endl;
}

void Route :: SetSelectionProbability(double sum){//Set selection probability for each route
  selection_probability = fitness/sum;
}

bool Route :: Check(){
  double check = (ncities-1)*ncities/2;
  double sum=0;
  for(int i=0; i<ncities; ++i){
    sum += genotype[i];
  }
  if (check == sum) return true;
  else return false;
}

void Route :: SetGenotype(std::vector<int> visiting_order){//Set order of visiting
  genotype = visiting_order;
}

void Route :: Permutation(int p1, int p2){//Permutation of two random city
  int temp;
  temp = genotype[p1];
  genotype[p1] = genotype[p2];
  genotype[p2] = temp;
}

void Route :: Shift(int places){//Shift of all cities
  rotate(genotype.rbegin(), genotype.rbegin()+places, genotype.rend());
}

void Route :: PartialShift(int start, int stop, int places){//Partial shift of a groupo of cities
  vector<int> shift;
  vector<int> old = genotype;
  rotate(genotype.begin(), genotype.begin()+start, genotype.end());
  for(int i=0; i<stop; ++i) shift.push_back(genotype[i]);
  genotype.erase(genotype.begin(), genotype.begin()+stop);
  genotype.insert(genotype.begin()+stop+places, shift.begin(), shift.end());
  rotate(genotype.rbegin(), genotype.rbegin()+start, genotype.rend());
  if (Check() == false ) genotype = old;
}

void Route :: Inversion(int start, int stop){//Inversion of a group of cities
  int temp;
  if (start>stop){
    temp = start;
    start = stop;
    stop = temp;
  }
  vector<int> inverted;
  vector<int> old = genotype;
  for(int i=start; i<stop; ++i) inverted.push_back(genotype[i]);
  reverse(inverted.begin(), inverted.end());
  genotype.insert(genotype.begin()+stop, inverted.begin(), inverted.end());
  genotype.erase(genotype.begin()+start, genotype.begin()+stop);
  if (Check() == false ) genotype = old;
}

void Route :: Print(int rank, int n){//Print route configuration to file
  ofstream Path;
  const int wd =12;
  int pos;
  Path.open("output/frames/path-"+to_string(rank)+"_"+to_string(n)+".xy");
  for(int i=0; i<ncities; ++i){
    pos = genotype[i];
    Path << setw(wd) << pos << setw(wd) << cities[pos].X() << setw(wd) << cities[pos].Y() << endl;
  }
  Path.close();
}

void Route :: Update(){//Update length of the route according to actual visiting order
  length = SetLength();
  fitness = Fitness();
}

bool operator<(const Route& a, const Route& b){//Sorting by decreasing fitness (increasing length)
  return a.Length() < b.Length();
}
