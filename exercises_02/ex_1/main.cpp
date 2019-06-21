/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019
Excercise 02.1

Francesco Ariele Piziali
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
//#include <string>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

double error(double sum, double sum2, int n);

int main (int argc, char *argv[]){

/*

  if (argc < 2) {
    cerr << "USAGE: " << argv[0] << " "<< endl;
  }

*/

/****************************************************************/
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

/****************************************************************/

  int M = 10000;
  int N = 100;
  int L = M/N;

  double r[M];
  double integral[N], integral2[N];
  double sum_prog[N], sum2_prog[N], err_prog[N];
  double sum;
  int k;

//Punto 1: uniform sampling
  ofstream output1("ex_1-1.out");

//Genero M numeri casuali tra 0 e 1
  for (int i=0; i<M; i++){
    r[i] = rnd.Rannyu();

  }



//Medio su N blocchi di lunghezza L=M/N
  for (int i=0; i<N; i++){
    sum = 0;
    for (int j=0; j<L; j++){
      k = j + i*L;
      sum += M_PI/2*cos(M_PI*0.5*r[k]);
    }
    integral[i] = sum/double(L);
    integral2[i] = integral[i]*integral[i];
  }

 //Somma cumulativa dei blocchi
  for (int i=0; i<N; i++){
    sum_prog[i] = 0;
    sum2_prog[i] = 0;
    err_prog[i] = 0;
  }

  for (int i=0; i<N; i++){
    for (int j=0; j<(i+1); j++){
      sum_prog[i] += integral[j];
      sum2_prog[i] += integral2[j];
    }
    sum_prog[i] /= double(i+1);
    sum2_prog[i] /= double(i+1);
    err_prog[i] = error(sum_prog[i], sum2_prog[i], i);
  }

  for (int i=0; i<N; i++){
    output1 << (i+1) << "\t" << sum_prog[i] << "\t" << err_prog[i] << endl;
  }
  output1.close();

//Punto 2: importance sampling
  ofstream output2("ex_1-2.out");

  double x, y;

  for (int i=0; i<M; i++){
    do{
      x = rnd.Rannyu();
      y = rnd.Rannyu();
    } while( y > 2*(1-x)/2) ;
      //(48/(24*M_PI-pow(M_PI,3))*(M_PI/2-pow(M_PI,3)/16*r[k]*r[k]))/(24*M_PI/(24*M_PI-pow(M_PI,3))));


    r[i] = x;

  }

  for (int i=0; i<N; i++){
    sum_prog[i] = 0;
    sum2_prog[i] = 0;
    err_prog[i] = 0;
    integral[i] = 0;
    integral2[i] = 0;
  }

  for (int i=0; i<N; i++){
    sum = 0;
    for (int j=0; j<L; j++){
      k = j + i*L;
      sum += (M_PI*0.5*cos(M_PI*0.5*r[k]))/(2*(1-r[k]));
      //r[k];
      //
      //(48/(24*M_PI-pow(M_PI,3))*(M_PI/2-pow(M_PI,3)/16*r[k]*r[k]));

      //(3/2*(1-r[k]*r[k]))*2/3;
    }
    integral[i] = sum/double(L);
    integral2[i] = integral[i]*integral[i];

  }

 //Somma cumulativa dei blocchi


  for (int i=0; i<N; i++){
    for (int j=0; j<(i+1); j++){
      sum_prog[i] += integral[j];
      sum2_prog[i] += integral2[j];
    }
    sum_prog[i] /= double(i+1);
    sum2_prog[i] /= double(i+1);
    err_prog[i] = error(sum_prog[i], sum2_prog[i], i);
  }

  for (int i=0; i<N; i++){
    output2 << (i+1) << "\t" << sum_prog[i] << "\t" << err_prog[i] << endl;
  }
  output2.close();

  rnd.SaveSeed();
  return 0;
}
/****************************************************************/

double error(double sum, double sum2, int n){
 if (n==0) return 0;
 else return (sqrt((sum2 - sum*sum)/n));
}



/****************************************************************
*****************************************************************
*****************************************************************
*****************************************************************/
