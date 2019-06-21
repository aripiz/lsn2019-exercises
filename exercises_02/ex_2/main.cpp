/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019
Excercise 02.2

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
double distance(double x, double y, double z);
double* blocking(int N_blocks, int L_elements, double *var, double *sum_prog, double *err_prog):


int main (int argc, char *argv[]){



  if (argc < 3) {
    cerr << "USAGE: " << argv[0] << " #Simulations #Steps"<< endl;
    return 1;
  }



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

  int M = atoi(argv[1]);

  int Nstep = atoi(argv[2]);

  double L_step = 1;

  double x[M], y[M], z[M];
  double r[Nstep][M];
  double r_max[Nstep], r_max2[Nstep], err[Nstep];
  int coordinate, direction;
  double sum,sum2;


//Punto 1: cubic lattice
  ofstream output1("ex_2-1.out");

  for (int i=0; i<Nstep; i++) {
    for (int j=0; j<M; j++) {
      r[i][j] = 0;
      x[j] = 0;
      y[j] = 0;
      z[j] = 0;
    }
  }

  r_max[0] = 0;
  err[0] = 0;
  for (int i=1; i<Nstep; i++){
    sum = 0;
    sum2 = 0;
    for (int j=0; j<M; j++){
      coordinate = rnd.RanInt(1,3);
      direction = rnd.RanInt(0,1);
      if (coordinate == 1){
        if (direction == 0) x[j] += L_step;
        else x[j] -= L_step;
      }
      if (coordinate == 2 ){
        if (direction == 0) y[j] += L_step;
        else y[j] -= L_step;
      }
      if (coordinate == 3){
        if (direction == 0) z[j] += L_step;
        else z[j] -= L_step;
      }
      r[i][j] = distance(x[j], y[j], z[j]);
      sum += r[i][j];
      sum2 += r[i][j]*r[i][j];
    }
    r_max[i] = sum/double(M);
    r_max2[i] = sum2/double(M);
    err[i] = error(r_max[i],r_max2[i], M);
  }


  output1 << M << "\t" << Nstep << endl;
  for (int i=0; i<Nstep; i++){
    output1 << (i+1) << "\t" << r_max[i] << "\t" << err[i] << endl;
  }

  output1.close();

  //Output RW singoli
  ofstream output3("cubicRWs.out");
  for (int i=0; i<Nstep; i++){
    output3 << i+1 << "\t";
    for (int j=0; j<10; j++)
      output3 << r[i][j] << "\t";
    output3 << endl;
  }
  output3.close();

//Punto 2: continuum
  ofstream output2("ex_2-2.out");
  double u, theta;

  for (int i=0; i<Nstep; i++) {
    for (int j=0; j<M; j++) {
      r[i][j] = 0;
      x[j] = 0;
      y[j] = 0;
      z[j] = 0;
    }
  }

  r_max[0] = 0;
  err[0] = 0;
  for (int i=1; i<Nstep; i++){
    sum = 0;
    sum2 = 0;
    for (int j=0; j<M; j++){
      //Random solid angle
      u = rnd.Rannyu(-1,1);
      theta = rnd.Rannyu(0,2*M_PI);
      x[j] += sqrt(1-u*u) * cos(theta);
      y[j] += sqrt(1-u*u) * sin(theta);
      z[j] += u;

      r[i][j] = distance(x[j], y[j], z[j]);
      sum += r[i][j];
      sum2 += r[i][j]*r[i][j];
    }
    r_max[i] = sum/double(M);
    r_max2[i] = sum2/double(M);
    err[i] = error(r_max[i], r_max2[i], M);
  }


  output2 << M << "\t" << Nstep << endl;
  for (int i=0; i<Nstep; i++){
    output2 << (i+1) << "\t" << r_max[i] << "\t" << err[i] << endl;
  }
  output2.close();


  //Output RW singoli
  ofstream output4("continuumRWs.out");
  for (int i=0; i<Nstep; i++){
    output4 << i+1 << "\t";
    for (int j=0; j<10; j++)
      output4 << r[i][j] << "\t";
    output4 << endl;
  }
  output4.close();


  rnd.SaveSeed();
  return 0;
}
/****************************************************************/
//FUNZIONI

double error(double sum, double sum2, int n){
 if (n==0) return 0;
 else return (sqrt((sum2 - sum*sum)/n));
}

double distance(double x, double y, double z){
  return sqrt(x*x + y*y + z*z);
}

double* blocking(int N_blocks, int L_elements, double *var, double *sum_prog, double *err_prog){
  double ave[N_blocks], ave2[N_blocks];
  double sum;
  double sum2_prog[N_blocks];
  int k;

  double final[2];

  for (int i=0; i<N_blocks; i++){
     sum = 0;
     for(int j=0; j<L_elements; j++){
       k = j + i*L;
       sum += var[k];
     }
     ave[i] = sum/double(L_elements);
     ave2[i] = (ave[i]*ave[i]);
  }

  for (int i=0; i<N_blocks; i++){
    for (int j=0; j<(i+1); j++){
      sum_prog[i] += ave[j];
      sum2_prog[i] += ave2[j];
    }
      sum_prog[i] /= (i+1);
      sum2_prog[i] /= (i+1);
      err_prog[i] = error(sum_prog[i], sum2_prog[i], i);
  }

  final[0] = sum_prog[N_blocks-1];
  final[1] = err_prog[N_blocks-1];
  return final;
}



/****************************************************************
*****************************************************************
*****************************************************************
*****************************************************************/
