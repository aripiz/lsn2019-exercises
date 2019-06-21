/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019
Excercise 08.1
Variational_Monte_Carlo.cpp
Francesco Ariele Piziali
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Variational_Monte_Carlo.h"

using namespace std;

int main()
{
  Input(); //Inizialization

  for(int ivar = 0; ivar <= nvar; ++ivar){
    if ((nvar>1) and (ivar>0)) Variation();

  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      Measure();
      Accumulate(); //Update block averages
      if((istep%10 == 0) && (print_xyz==true)){
        ConfXYZ();//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
      }
    }
    Averages(iblk);   //Print results for current block
  }
  if (nvar>1) Final(ivar);
}
  return 0;
}


void Input(void)
{
  ifstream ReadInput,ReadConf;
  cout << "Variational Monte Carlo simulation             " << endl << endl;
  cout << "Ground state of 1D quantum particle in the external potential:" << endl;
  cout << "   V(x) = x^4 - 2.5x^2" << endl << endl;
  cout << "Test wave function:" << endl;
  cout << "   Psi_T(x) = exp[ -(x-mean)^2/2sigma^2 ] + exp[ -(x+mean)^2/2sigma^2 ]" << endl << endl;


//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> x;

  ReadInput >> mean;

  ReadInput >> sigma;

  ReadInput >> delta;

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> delta_var;

  ReadInput >> nvar;


  ReadInput >> print_xyz;

  ReadInput >> print_instant;

  ReadInput >> print_blocks;

  cout << "Parameters of Psi_T(x): " << endl;
  cout << "   mean = " << mean << endl << "   sigma = " << sigma << endl << endl;

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

//Prepare arrays for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic
  ie = 2; //Total

  n_props = 3; //Number of observables

//Initial configuration

//Evaluate potential energy and virial of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Starting position = " << x << endl;
  cout << "Initial potential energy  = " << walker[iv] << endl << endl;
  accepted_var = 0.0;
  attempted_var = 0.0;

}

void Move(void)
{
  double p, p_old, p_new;
  double xold,xnew;

  //Old
    xold = x;
    p_old = PsiT(xold);

  //New
    xnew = x + delta*rnd.Rannyu(-1,1);
    p_new = PsiT(xnew);

  //Metropolis test
    p = p_new/p_old;
    if(p >= rnd.Rannyu())
    {
    //Update
       x = xnew;

       accepted = accepted + 1.0;
    }
    attempted = attempted + 1.0;
}

double PsiT(double x){

  return exp( -(x-mean)*(x-mean)/(sigma*sigma) ) + exp( -(x+mean)*(x+mean)/(sigma*sigma) ) + 2*exp( -(x*x+mean*mean)/(sigma*sigma) );

}

double Potential(double x){
  return x*x*x*x - 2.5*x*x;
}

double Kinetic(double x){
  double m = (x-mean)*(x-mean)/(sigma*sigma);
  double p = (x+mean)*(x+mean)/(sigma*sigma);
  return -0.5 * ( (m-1)*exp(-m/2)  + (p-1)*exp(-p/2) ) / (sigma*sigma);
}

void Measure()
{
  double v, k;

// contribution to energies
  v = Potential(x);
  k = Kinetic(x);

//update walker
  walker[iv] = v;
  walker[ik] = k;
  walker[ie] = k+v;

  //Print instant values
  if (print_instant == true) {
    ofstream InstantEkin, InstantEpot;

    InstantEpot.open("output/instant.epot.0",ios::app);
    InstantEkin.open("output/instant.ekin.0",ios::app);

    InstantEpot << walker[iv] << endl; //Potential energy
    InstantEkin << walker[ik] << endl; //Kinetic energy

    InstantEpot.close();
    InstantEkin.close();
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
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) //Print results for current block
{
   ofstream Epot, Ekin, Etot;
   const int wd=12;


    stima[iv] = blk_av[iv]/blk_norm; //Potential energy
    glob_av[iv] += stima[iv];
    glob_av2[iv] += stima[iv]*stima[iv];
    err[iv]=Error(glob_av[iv],glob_av2[iv],iblk);

    stima[ik]  = blk_av[ik]/blk_norm; //Potential energy
    glob_av[ik] += stima[ik];
    glob_av2[ik] += stima[ik]*stima[ik];
    err[ik]=Error(glob_av[ik],glob_av2[ik],iblk);

    stima[ie] = blk_av[ie]/blk_norm; //Potential energy
    glob_av[ie] += stima[ie];
    glob_av2[ie] += stima[ie]*stima[ie];
    err[ie]=Error(glob_av[ie],glob_av2[ie],iblk);


    if (print_blocks == true) {
      Epot.open("output/output.epot.0", ios::app);
      Ekin.open("output/output.ekin.0", ios::app);
      Etot.open("output/output.etot.0", ios::app);

//Potential energy per particle
      Epot << setw(wd) << iblk <<  setw(wd) << stima[iv] << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err[iv] << endl;
//Kinetic
      Ekin << setw(wd) << iblk <<  setw(wd) << stima[ik] << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err[ik] << endl;
//Total
      Etot << setw(wd) << iblk <<  setw(wd) << stima[ie] << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err[ie] << endl;


      cout << "Block number " << iblk << endl;
      cout << "Acceptance rate " << accepted/attempted << endl << endl;
      cout << "----------------------------" << endl << endl;

      Epot.close();
      Ekin.close();
      Etot.close();
    }

}

void ConfXYZ(){ //Write configuration in .xyz format
  ofstream WriteX;

  WriteX.open("output/config.final", ios::app);
  WriteX << x << endl;

  WriteX.close();
}

void Final(int ivar)
{
  const int wd = 12;
  energy_new = glob_av[ie]/(double)nblk;

  if (ivar >= 1){
    cout << "Step number " << ivar << endl;
    cout << "Acceptance rate " << accepted_var/attempted_var << endl;
    cout << "Delta " << delta_var << endl << endl;

    cout << "Final values: " << endl;
    cout << setw(wd) << "mean = " << setw(wd) << mean << endl;
    cout << setw(wd) << "sigma = " << setw(wd) << sigma << endl;
    cout << setw(wd) << "<H> = " << setw(wd) << energy_new  <<" +/- " << err[ie] << endl;
  }
  ofstream Var;

  Var.open("output/output.var.0", ios::app);
  Var << setw(wd) << accepted_var << setw(wd) << mean <<  setw(wd) << sigma << setw(wd) << energy_new << setw(wd) << err[ie] << endl;
  Var.close();

  if (ivar == 0){
    energy_old = energy_new;
    energy_opt = energy_new;
    mean_opt = mean;
    sigma_opt = sigma;
    err_opt = err[ie];
  }

  if ((energy_new < energy_old)) {
    cout << "<H> decreased: new parameters accepted!" << endl;
    energy_old = energy_new;
    delta_var /= 1.25;
    accepted_var += 1.0;

    energy_opt = energy_new;
    mean_opt = mean;
    sigma_opt = sigma;
    err_opt = err[ie];
  }
    if (energy_new > energy_old){
    cout << "<H> increased: new parameters refused!" << endl;
    mean = mean_old;
    sigma = sigma_old;
  }
  cout << endl;
  cout << "Optimal values: " << endl;
  cout << setw(wd) << "mean = " << setw(wd) << mean_opt << endl;
  cout << setw(wd) <<"sigma = " << setw(wd) << sigma_opt << endl;
  cout << setw(wd) << "<H> = " << setw(wd) << energy_opt << " +/- " << err_opt << endl;

  cout << "----------------------------" << endl << endl;

}

void Variation(){

  mean_old = mean;
  sigma_old = sigma;
  //int sign = 1;
  //if (rnd.Rannyu() < 0.5) sign = -1;
  if (rnd.Rannyu() < 0.5) mean += delta_var * rnd.Rannyu(-1,1);//*sign;
  else sigma += delta_var * rnd.Rannyu(-1,1);//*sign;

  attempted_var += 1.0;

}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

/****************************************************************
*****************************************************************
*****************************************************************
*****************************************************************/
