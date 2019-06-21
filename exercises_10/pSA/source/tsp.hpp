/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019

Travelling Salesman Problem
tsp.hpp

Francesco Ariele Piziali
*****************************************************************
*****************************************************************/

#include<vector>

#ifndef __Tsp__
#define __Tsp__

class City {
private:
  double x,y;

public:
  City(){};
  City(double a, double b) : x(a), y(b) {};
  ~City(){};
  double X(){return x;}
  double Y(){return y;}
  double Distance(City);
  void Show();

};



class Route {
private:
  int ncities;
  std::vector<City> cities;
  std::vector<int> genotype;
  double length, fitness, selection_probability;

protected:
  double SetLength();
  int Pbc(int);

public:
  Route(){};
  Route(std::vector<City>, std::vector<int>);
  ~Route(){};

  void Show();
  double Fitness() const {return 1/length;}
  double Length() const {return length;}
  void SetSelectionProbability(double);
  double SelectionProbability(){return selection_probability;}
  void Exclude(){selection_probability = 0;}
  std::vector<int> Genotype(){return genotype;}
  bool Check();

  void Print(int, int);
  void Update();
  void SetGenotype(std::vector<int>);

  void Permutation(int, int);
  void Shift(int);
  void PartialShift(int,int,int);
  void Inversion(int, int);
};

bool operator<(const Route& a, const Route& b);


#endif
