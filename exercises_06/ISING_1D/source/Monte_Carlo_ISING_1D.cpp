/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main(int argc, char *argv[])
{
  if(argc>1) input_temp = atof(argv[1]);
  Input(); //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();

//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  if (input_temp>0) temp = input_temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;

  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> read_initial;
  ReadInput >> instant_print;
  ReadInput >> verbose;

  if(metro==1) cout << "The program performs Metropolis moves" << endl;
  else cout << "The program performs Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility

  n_props = 4; //Number of observables

//Prepare folder
  if (metro == true) folder = string("./metropolis_") + to_string(h)+string("/");
  else folder = string("./gibbs_") + to_string(h) +string("/");
  //cout << folder + string("averages.out")<< endl;

  //read initial configuration from file
  if (read_initial == true ){
    ifstream InputConf;
    InputConf.open("config.final");
    if (InputConf.is_open()){
      for (int i=0; i<nspin; ++i) InputConf >> s[i];
      InputConf.close();
      cout << "Initial configuration loaded from file config.final" << endl;
    }
    else{
      cout << "File config.final not found: preparing random initial configuration " << endl;
      read_initial = false;
    }
  }

  //random initial configuration
  if (read_initial == false  ){
    for (int i=0; i<nspin; ++i)
    {
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
    }
    cout << "Preparing random initial configuration " << endl;
  }


//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}

void Move(int metro)
{
  int o;
  double s_old, energy_old, energy_new, s_new;
  //double energy_up, energy_down;
  double alpha;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==true) //Metropolis
    {
      ++attempted;
      s_old = s[o];
      energy_old = Boltzmann(s_old,o);
      if (s_old == +1) s_new = -1;
      else s_new = +1;
      energy_new = Boltzmann(s_new,o);
      alpha = exp(-beta*(energy_new-energy_old));
      if ( rnd.Rannyu()<fmin(1, alpha) ) {
        s[o] = s_new;
        ++accepted;
      }
    }

    else //Gibbs sampling
    {
      ++attempted;
      alpha = 1/(1+exp(-2*beta*J*(s[Pbc(o-1)] + s[Pbc(o+1)])-2*beta*h));
      if ( rnd.Rannyu()<alpha ) s_new = +1;
      else s_new = -1;

      s[o] = s_new;
      ++accepted;

    }
  }
}

double Boltzmann(double sp, int ip)
{
  double ene = -J * sp * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sp;
  return ene;
}

void Measure()
{
  //int bin;
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);

     m += s[i];

  }
  walker[iu] = u;

// INCLUDE YOUR CODE HERE
  walker[im] = m;
  walker[ic] = u*u;
  walker[ix] = m*m;

  if (instant_print == true) {
    ofstream InstantMagn, InstantEpot;

    InstantEpot.open("output/instant.epot.0",ios::app);
    InstantMagn.open("output/instant.magn.0",ios::app);

    InstantEpot << walker[iu]/(double)nspin << endl; //Instant potential energy
    InstantMagn << walker[im]/(double)nspin << endl; //Instant Magnetization

    InstantEpot.close();
    InstantMagn.close();
  }


}

void Reset(int iblk) //Reset block averages
{

   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}

void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] += walker[i];
   }
   blk_norm += 1.0;
}

void Averages(int iblk) //Print results for current block
{

   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;

    Ene.open("output/output.ene.0",ios::app);

    stima[iu] = blk_av[iu]/blk_norm/(double)nspin; //Energy

    glob_av[iu]  += stima[iu];
    glob_av2[iu] += stima[iu]*stima[iu];
    err[iu]=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima[iu] << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err[iu] << endl;
    Ene.close();

// INCLUDE YOUR CODE HERE
    Mag.open("output/output.mag.0",ios::app);

    stima[im] = blk_av[im]/blk_norm/(double)nspin; //Magnetization

    glob_av[im]  += stima[im];
    glob_av2[im] += stima[im]*stima[im];
    err[im]=Error(glob_av[im],glob_av2[im],iblk);
    Mag << setw(wd) << iblk <<  setw(wd) << stima[im] << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err[im] << endl;
    Mag.close();

    Heat.open("output/output.heat.0",ios::app);

    stima[ic] = beta*beta * (blk_av[ic]/blk_norm - blk_av[iu]/blk_norm*blk_av[iu]/blk_norm)/(double)nspin; //Heat capacity

    glob_av[ic]  += stima[ic];
    glob_av2[ic] += stima[ic]*stima[ic];
    err[ic]=Error(glob_av[ic],glob_av2[ic],iblk);
    Heat << setw(wd) << iblk <<  setw(wd) << stima[ic] << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err[ic] << endl;
    Heat.close();

    Chi.open("output/output.chi.0",ios::app);

    stima[ix] = beta * (blk_av[ix]/blk_norm - stima[im]*stima[im])/(double)nspin; //Susceptibility

    glob_av[ix]  += stima[ix];
    glob_av2[ix] += stima[ix]*stima[ix];
    err[ix]=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi << setw(wd) << iblk <<  setw(wd) << stima[ix] << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err[ix] << endl;
    Chi.close();

    if (verbose ==true ){
     cout << "Block number " << iblk << endl;
     cout << "Acceptance rate " << accepted/attempted << endl << endl;
     cout << "----------------------------" << endl << endl;
   }


}


void ConfFinal(void)
{
  ofstream WriteConf, FinalAverages;
  const int wd=12;
  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  cout << "Print final averages to file averages.out" << endl << endl;
  FinalAverages.open("output/averages.out", ios::app);
  FinalAverages << setw(wd) << temp ;
  FinalAverages << setw(wd) << h;
  FinalAverages << setw(wd) << glob_av[iu]/(double)nblk << setw(wd) << err[iu];
  FinalAverages << setw(wd) << glob_av[im]/(double)nblk << setw(wd) << err[im];
  FinalAverages << setw(wd) << glob_av[ic]/(double)nblk << setw(wd) << err[ic];
  FinalAverages << setw(wd) << glob_av[ix]/(double)nblk << setw(wd) << err[ix];
  FinalAverages << endl;
  FinalAverages.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk -pow(sum/(double)iblk,2))/(double)iblk);
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
