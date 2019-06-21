/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __fluid_
#define __fluid_

//Random numbers
#include "random.h"
int seed[4];
Random rnd("seed.in");

//parameters, observables
const int m_props=1000;
int n_props, iv, ik, ie;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props], err[m_props];
double stima[m_props];

//configuration
double x;

//parameters
double mean, sigma;

// simulation
int nstep, nblk;
double delta;
bool print_xyz, print_instant, print_blocks;

//variation
double energy_new, energy_old, energy_opt, err_opt;
int nvar;
double delta_var, accepted_var, attempted_var;
double mean_old, sigma_old, mean_opt, sigma_opt;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfXYZ();
void Measure(void);
double Potential(double);
double PsiT(double);
void Final(int);
double Error(double,double,int);
void Variation();

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
