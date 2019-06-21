/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019
Excercise 03.1

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
double* blocking(int N_blocks, int L_elements, double *var, double *sum_prog, double *err_prog);



int main (int argc, char *argv[]){

//Parametri in ingresso
  if (argc < 6) {
    cerr << "USAGE: " << argv[0] << " [asset price] [delivery time] [strike price] [risk-free rate] [volatility]"<< endl;
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

  //Numero di simulazioni e blocchi
  int M = 10000;
  int N = 100;
  int L = M/N;

  //GBM parameters
  double r = atof(argv[4]);
  double sigma = atof(argv[5]);
  double S_0 = atof(argv[1]);
  double K = atof(argv[3]);
  double T = atof(argv[2]);

  //Asset prices, options
  double S_t[M];
  double W;
  double C[M], P[M];

  //Variabili per blocking stats
  double C_ave[N], P_ave[N], C_ave2[N], P_ave2[N];
  double sum_C, sum_P;
  double err_prog_C[N], err_prog_P[N];
  double sum_prog_C[N], sum_prog_P[N];
  double sum2_prog_C[N], sum2_prog_P[N];
  int k;


  //Punto 1: Direct
  for (int i=0; i<M; i++){
    W = rnd.Gauss(0,T);
    S_t[i] = S_0 * exp((r-0.5*sigma*sigma)*T+sigma*W);
    C[i] = exp(-r*T) * fmax(0, S_t[i] - K);
    P[i] = exp(-r*T) * fmax(0, K - S_t[i]);
  }

  //Blocking stat
  for (int i=0; i<N; i++){
     sum_C = 0;
     sum_P = 0;
     for(int j=0; j<L; j++){
       k = j + i*L;
       sum_C += C[k];
       sum_P += P[k];
     }
     C_ave[i] = sum_C/double(L);
     C_ave2[i] = (C_ave[i]*C_ave[i]);
     P_ave[i] = sum_P/double(L);
     P_ave2[i] = (P_ave[i]*P_ave[i]);
  }

  for (int i=0; i<N; i++){
    for (int j=0; j<(i+1); j++){
      sum_prog_C[i] += C_ave[j];
      sum2_prog_C[i] += C_ave2[j];

      sum_prog_P[i] += P_ave[j];
      sum2_prog_P[i] += P_ave2[j];
    }
      sum_prog_C[i] /= (i+1);
      sum2_prog_C[i] /= (i+1);
      err_prog_C[i] = error(sum_prog_C[i], sum2_prog_C[i], i);

      sum_prog_P[i] /= (i+1);
      sum2_prog_P[i] /= (i+1);
      err_prog_P[i] = error(sum_prog_P[i], sum2_prog_P[i], i);
  }

  ofstream output_C("call.out");
  ofstream output_P("put.out");

  for (int i=0; i<N; i++){
    output_C << (i+1) << "\t" << sum_prog_C[i] << "\t" << err_prog_C[i] << endl;
    output_P << (i+1) << "\t" << sum_prog_P[i] << "\t" << err_prog_P[i] << endl;
  }
  output_C.close();
  output_P.close();

  //Punto 2: Discrete
  for (int i=0; i<N; i++){
    C_ave[i] = 0;
    C_ave2[i] = 0;
    P_ave[i] = 0;
    P_ave2[i] = 0;
    sum_prog_C[i] = 0;
    sum2_prog_C[i] = 0;
    sum_prog_P[i] = 0;
    sum2_prog_P[i] = 0;
  }

  double S_discrete;
  int steps = 100;
  double delta_t = T/steps;

  for (int i=0; i<M; i++){
    S_discrete = S_0;
    for (int j=0; j<steps; j++){
      W = rnd.Gauss(0,T);
      S_discrete *= exp((r-0.5*sigma*sigma)*delta_t+sigma*W*sqrt(delta_t));
    }
    S_t[i] = S_discrete;
    C[i] = exp(-r*T) * fmax(0, S_t[i] - K);
    P[i] = exp(-r*T) * fmax(0, K - S_t[i]);
    //cout << P[i] << endl;
  }

  //Blocking stat
  for (int i=0; i<N; i++){
     sum_C = 0;
     sum_P = 0;
     for(int j=0; j<L; j++){
       k = j + i*L;
       sum_C += C[k];
       sum_P += P[k];
     }
     C_ave[i] = sum_C/double(L);
     C_ave2[i] = (C_ave[i]*C_ave[i]);
     P_ave[i] = sum_P/double(L);
     P_ave2[i] = (P_ave[i]*P_ave[i]);
  }

  for (int i=0; i<N; i++){
    for (int j=0; j<(i+1); j++){
      sum_prog_C[i] += C_ave[j];
      sum2_prog_C[i] += C_ave2[j];

      sum_prog_P[i] += P_ave[j];
      sum2_prog_P[i] += P_ave2[j];
    }
      sum_prog_C[i] /= (i+1);
      sum2_prog_C[i] /= (i+1);
      err_prog_C[i] = error(sum_prog_C[i], sum2_prog_C[i], i);

      sum_prog_P[i] /= (i+1);
      sum2_prog_P[i] /= (i+1);
      err_prog_P[i] = error(sum_prog_P[i], sum2_prog_P[i], i);
  }

  ofstream output_C2("call_discrete.out");
  ofstream output_P2("put_discrete.out");
  for (int i=0; i<N; i++){
    output_C2 << (i+1) << "\t" << sum_prog_C[i] << "\t" << err_prog_C[i] << endl;
    output_P2 << (i+1) << "\t" << sum_prog_P[i] << "\t" << err_prog_P[i] << endl;
  }
  output_C2.close();
  output_P2.close();

  //Output parametri usati
  ofstream output_par("parameters.out");
  output_par << "simulations" << "\t" << M << endl;
  output_par << "asset price" << "\t" << S_0 << endl;
  output_par << "delivery time" << "\t" << T << endl;
  output_par << "strike price" <<  "\t" << K << endl;
  output_par << "risk-free rate" << "\t" << r << endl;
  output_par << "volatility" << "\t" << sigma << endl;
  output_par.close();

  rnd.SaveSeed();
  return 0;
}
/****************************************************************/
//FUNZIONI

double error(double sum, double sum2, int n){
 if (n==0) return 0;
 else return (sqrt((sum2 - sum*sum)/n));
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
