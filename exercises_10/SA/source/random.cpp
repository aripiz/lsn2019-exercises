/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019
Random.cpp

Francesco Ariele Piziali
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"


using namespace std;

/*****************************************************************/
// Methods of Random

Random :: Random(std::string seed_file){
  prime = 0;
  ifstream Primes("Primes");
  if (Primes.is_open()){
    Primes >> p1 >> p2 ;
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();
  ifstream input(seed_file);
  std::string property;
  if ( input.is_open() ){
    while ( !input.eof() ){
       input >> property;
       if( property == "RANDOMSEED" ){
          input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
          SetRandom(seed,p1,p2);
       }
    }
    input.close();
  } else cerr << "PROBLEM: Unable to open " << seed_file << endl;
}

Random :: Random(std::string seed_file, int set){
  prime = set;
  ifstream Primes("Primes");
  int count = 0;
  if (Primes.is_open()){
    do{
      Primes >> p1 >> p2 ;
      count++;
    } while(count>=prime);
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();
  ifstream input(seed_file);
  std::string property;
  if ( input.is_open() ){
    while ( !input.eof() ){
       input >> property;
       if( property == "RANDOMSEED" ){
          input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
          SetRandom(seed,p1,p2);
       }
    }
    input.close();
  } else cerr << "PROBLEM: Unable to open " << seed_file << endl;
}


Random :: ~Random(){
  SaveSeed();
}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed_"+ to_string(prime)+ ".out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open seed.out" << endl;
  WriteSeed.close();
  return;
}

int Random :: Integer(int min, int max){
  double y = Rannyu(min, max+1);
  while ( y == (max+1) ) y = Rannyu(min, max+1);
  return int(y);
}

double Random :: Expo(double rate) {
  double y = Rannyu();
  while ( y == 1 ) y = Rannyu();
  return (-1/rate * log(1-y));
}

double Random :: CauchyLorentz(double mean, double gamma){
  double y = Rannyu();
  return mean + gamma * tan(M_PI*(y-1/2));
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0]%4096;
  l2 = s[1]%4096;
  l3 = s[2]%4096;
  l4 = s[3]%4096;
  l4 = 2*(l4/2)+1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}
/****************************************************************/






/****************************************************************
*****************************************************************
*****************************************************************
*****************************************************************/
