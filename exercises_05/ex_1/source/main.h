/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019
Excercise 05.1
main.h
Francesco Ariele Piziali
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include "random.h"

using namespace std;

//Simulation parameters
int M; //Number of steps
const int N = 100; //Number of blocks
double delta; //Step lenght

//Coordinates
double x_old, y_old, z_old; //Old coordinates
double x_new, y_new, z_new; //New coordinates

vector<double> r;  // Distances from the origin of the sampled points

double p_old, p_new; //p(r)

//Random moves
double u, theta; //Random solid angle
double tx,ty,tz; //Normal transition

//Acceptance rate
int accepted; //Number of accep moves
bool accept; //Result of acceptance test

//Output
ofstream CoorOut;
ofstream AccOut;

/*****************************************************************/
//Functions
double error(double sum, double sum2, int n);
double distance(double x, double y, double z);
void blocking(vector<double> var, int nblocks, string output_file);
bool test(double p_new, double p_old, double t);

//Hydrogen orbitals probability functions
double p_1s(double x, double y, double z);
double p_2p(double x, double y, double z);
double p_320(double x, double y, double z);
double p_321(double x, double y, double z);
