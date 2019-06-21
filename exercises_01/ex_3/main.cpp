/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019
Excercise 01.3

Francesco Ariele Piziali
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

double error(double sum, double sum2, int n){
  if (n==0) return 0;
  else return (sqrt((sum2 - sum*sum)/n));
}


int main (int argc, char *argv[]){

  Random rnd;
  int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()){
    Primes >> p1 >> p2 ;
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();
  ifstream input("seed.in");
  string property;
  if (input.is_open()){
    while ( !input.eof() ){
       input >> property;
       if( property == "RANDOMSEED" ){
          input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
          rnd.SetRandom(seed,p1,p2);
       }
    }
    input.close();
  } else cerr << "PROBLEM: Unable to open seed.in" << endl;

/**************************************************************/

  int M = 100000;
  int N = 100;
  int L = M/N;

  double d = 1.0;
  double l = 0.25 * d;

  int N_hit;
  int k;
  double y_center[M], phi[M];
  double x_angle, y_angle;
  double pi[N], pi2[N];
  double sum_prog[N], sum2_prog[N], err_prog[N];

  ofstream output1("ex_3.out");
  /**************************************************************/

  //Campiono ordinata del centro dell'ago [0,d/2] e
  //angolo rispetto all'orizzontale [0,pi/2]
  for (int i=0; i<M; i++){
    y_center[i] = rnd.Rannyu(0, d/2);
    //Per evitare di usare M_PI
    do{
      x_angle = rnd.Rannyu();
      y_angle= rnd.Rannyu();
    } while (x_angle*x_angle + y_angle*y_angle > 1);
    phi[i] = atan(y_angle/x_angle);
    //phi[i] = rnd.Rannyu(0, M_PI/2);
  }

  //Conto quanti aghi attraversano le linee e calcolo pi
  for (int i=0; i<N; i++){
    N_hit = 0;
    for (int j=0; j<L; j++){
      k = j + i*L;
      if ( (l/2*sin(phi[k])) >= (y_center[k])) N_hit++;
    }
    pi[i] =(2*l*L)/(N_hit*d);
    pi2[i] = (pi[i]*pi[i]);
    //sum += pi[i];
  }

  //Statistica pi
  for (int i=0; i<N; i++){
    sum_prog[i] = 0;
    sum2_prog[i] = 0;
  }

  for (int i=0; i<N; i++){
    for (int j=0; j<(i+1); j++){
      sum_prog[i] += pi[j];
      sum2_prog[i] += pi2[j];
    }
    sum_prog[i] /= (i+1);
    sum2_prog[i] /= (i+1);
    err_prog[i] = error(sum_prog[i], sum2_prog[i], i);
  }

  for (int i=0; i<N; i++){
    output1 << (i+1) << "\t" << sum_prog[i] << "\t" << err_prog[i] << endl;
  }

  output1.close();

  rnd.SaveSeed();
  return 0;
}

/****************************************************************
*****************************************************************
*****************************************************************
*****************************************************************/
