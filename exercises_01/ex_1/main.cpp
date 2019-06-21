/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019
Excercise 01.1

Francesco Ariele Piziali
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

double error(double sum, double sum2, int n);


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



  int M = 1000000;
  int N = 100;
  int L = M/N;

  double r[M];
  double ave[N], ave2[N];
  double sum_prog[N], sum2_prog[N], err_prog[N];
  double sum;
  int k;

//Punto 1: media e incertezza
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
      sum += r[k];
    }
    ave[i] = sum/double(L);
    ave2[i] = (ave[i]*ave[i]);
   }

   //Somma cumulativa dei blocchi
   for (int i=0; i<N; i++){
     for (int j=0; j<(i+1); j++){
       sum_prog[i] += ave[j];
       sum2_prog[i] += ave2[j];
     }
       sum_prog[i] /= (i+1);
       sum2_prog[i] /= (i+1);
       err_prog[i] = error(sum_prog[i], sum2_prog[i], i);
   }

   for (int i=0; i<N; i++){
     output1 << (i+1) << ";" << sum_prog[i] << ";" << err_prog[i] << endl;
   }

   output1.close();


//Punto 2: dev. std. e incertezza

   for (int i=0; i<N; i++){
     ave[i] = 0;
     ave2[i] = 0;
     sum_prog[i] = 0;
     sum2_prog[i] = 0;
     err_prog[i] = 0;
   }

   ofstream output2("ex_1-2.out");

   for (int i=0; i<N; i++){
      sum = 0;
      for(int j=0; j<L; j++){
        k = j + i*L;
        sum += (r[k] - 0.5)*(r[k] - 0.5);
      }
      ave[i] = sum/double(L);
      ave2[i] = (ave[i]*ave[i]);
   }

   for (int i=0; i<N; i++){
     for (int j=0; j<(i+1); j++){
       sum_prog[i] += ave[j];
       sum2_prog[i] += ave2[j];
     }
       sum_prog[i] /= (i+1);
       sum2_prog[i] /= (i+1);
       err_prog[i] = error(sum_prog[i], sum2_prog[i], i);
   }

   for (int i=0; i<N; i++){
     output2 << (i+1) << ";" << sum_prog[i] << ";" << err_prog[i] << endl;
   }

   output2.close();

//Punto 3: chi-quadro
  ofstream output3("ex_1-3.out");

  //Divido l'intervallo [0,1] in N_bin sottointervalli
  int N_bin = 100;
  int n = M/N;
  int bin[N_bin];
  double chiq[N],chiq2[N];
  for (int i=0; i<N; i++) {
    chiq[i] = 0;
    chiq2[i] = 0;
  }

  //Calcolo k-esimo chi-quadro
  for(int k=0; k<N; k++){
    for (int i=0; i<N_bin; i++) bin[i] = 0; //Azzero contatori bin
    //Metto dati nei bin
    for (int i=k*n; i<(k*n+n); i++){
      for (int j=0; j<N_bin; j++){
        if ( (r[i] >= j/double(N_bin)) && (r[i] < (j+1)/double(N_bin)) ) {
          bin[j]++;
          break;
        }
      }
    }
    for (int i=0; i<N; i++){
      chiq[k] += (bin[i] - n/double(N_bin))*(bin[i] - n/double(N_bin)) / (n/double(N_bin));
      chiq2[k] = chiq[k]*chiq[k];
    }
  }

  //Statistica chi-quadro
  for (int i=0; i<N; i++){
    sum_prog[i] = 0;
    sum2_prog[i] = 0;
    err_prog[i] = 0;
  }

  for (int i=0; i<N; i++){
    for (int j=0; j<(i+1); j++){
      sum_prog[i] += chiq[j];
      sum2_prog[i] += chiq2[j];
    }
      sum_prog[i] /= (i+1);
      sum2_prog[i] /= (i+1);
      err_prog[i] = error(sum_prog[i], sum2_prog[i], i);
  }


  for(int i=0; i<N; i++){
    output3 << i+1 << ";" << sum_prog[i] << ";" << err_prog[i] << endl;
  }

  output3.close();


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
