/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019
Excercise 01.2

Francesco Ariele Piziali
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;


int main (int argc, char *argv[]){

/*  */

  if (argc < 2) {
    cerr << "USAGE: " << argv[0] << " Number of random variable to sum"<< endl;
  }


/*  */

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

  string output_file;

  int N = atoi(argv[1]);
  int L =  10000;
  double sum;
  double ave[L];
  //double sum_prog[N], sum2_prog[N], err_prog[N];
  //int k;

/*
  //Genero numeri casuali secondo le 3 distr
  int dice[M];
  double exp[M];
  double lor[M];
  for (int i=0; i<M; i++){
    dice[i] = rnd.RanInt(1,6);
    exp[i] = rnd.Expo(1);
    lor[i] = rnd.CauchyLorentz(0,1);
  }
*/




	//Figura1: Dado standard
  for (int i=0; i<L; i++) ave[i] = 0;
  for (int i=0; i<L; i++){
    sum = 0;
    for (int j=0; j<N; j++){
      sum += rnd.RanInt(1,6);
    }
    ave[i] = sum/double(N);
  }

  output_file = string("std_dice_N-") + argv[1] + string(".out") ;
  ofstream output1(output_file);
  for (int i=0; i<L; i++){
    output1 << i+1 << "\t" << ave[i] << endl;
  }
  output1.close();

  //Figura2: Esponenziale
  for (int i=0; i<L; i++) ave[i] = 0;
  for (int i=0; i<L; i++){
    sum = 0;
    for (int j=0; j<N; j++){
      sum += rnd.Expo(1);
    }
    ave[i] = sum/double(N);
  }

  output_file = string("exp_dice_N-") + argv[1] + string(".out") ;
	ofstream output2(output_file);
  for (int i=0; i<L; i++){
    output2 << i+1 << "\t" << ave[i] << endl;
  }
  output2.close();

  //Figura3: Lorentziana
  for (int i=0; i<L; i++) ave[i] = 0;
  for (int i=0; i<L; i++){
    sum = 0;
    for (int j=0; j<N; j++){
      sum += rnd.CauchyLorentz(0,1);
    }
    ave[i] = sum/double(N);
  }

  output_file = string("lor_dice_N-") + argv[1] + string(".out") ;
  ofstream output3(output_file);
  for (int i=0; i<L; i++){
    output3 << i+1 << "\t" << ave[i] << endl;
  }
  output3.close();



  rnd.SaveSeed();
  return 0;
}

/****************************************************************
*****************************************************************
*****************************************************************
*****************************************************************/
