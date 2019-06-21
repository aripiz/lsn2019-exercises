/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019
Excercise 04.1
Molecular_Dynamics_NVE.h
Francesco Ariele Piziali
*****************************************************************
*****************************************************************/

#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>

#include "Molecular_Dynamics_NVE.h"
#include "random.h"


using namespace std;

int main(){
  Input(); //Inizialization
  int nconf = 1;
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
      if(istep%10 == 0){
         Measure();     //Properties measurement
         if(print_xyz == true) ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
         nconf += 1;
         Accumulate(); //Update block averages
      }
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration
  return 0;
}

//Functions
void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf,ReadOld;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

   //Initialize random number generator
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

  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;

  ReadInput >> nstep;
  ReadInput >> nblk;

  ReadInput >> iprint;
  ReadInput >> old;
  ReadInput >> rescale;
  ReadInput >> print_xyz;
  ReadInput >> instant_print;


  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Tail corrections for potential energy and pressure
  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
  vtail = (8.0*M_PI*rho)/(9.0*pow(rcut,9)) - (8.0*M_PI*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*M_PI*rho)/(9.0*pow(rcut,9)) - (16.0*M_PI*rho)/(3.0*pow(rcut,3));
  cout << "Tail correction for the potential energy = " << vtail << endl;
  cout << "Tail correction for the virial           = " << ptail << endl<<endl;

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  iw = 4, //Virial
  n_props = 5; //Number of observables

  //measurement of g(r)
    igofr = n_props;
    nbins = 100;
    n_props = n_props + nbins;
    bin_size = (box/2.0)/(double)nbins;

  //Resize vectors based on npart and n_props
  x.resize(npart);
  y.resize(npart);
  z.resize(npart);
  xold.resize(npart);
  yold.resize(npart);
  zold.resize(npart);
  vx.resize(npart);
  vy.resize(npart);
  vz.resize(npart);

  walker.resize(n_props);
  blk_av.resize(n_props);
  glob_av.resize(n_props);
  glob_av2.resize(n_props);
  err.resize(n_props);

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

  //Read old configuration (if there are)
  if (old == true){
    ReadOld.open("old.0");

    if (ReadOld.is_open()) {
      cout << "Read old configuration from file old.0 " << endl << endl;
      for (int i=0; i<npart; ++i){
        ReadOld >> xold[i] >> yold[i] >> zold[i];
        xold[i] = xold[i] * box;
        yold[i] = yold[i] * box;
        zold[i] = zold[i] * box;
      }
      ReadOld.close();

      //Rescaling
      if (rescale == true){
        cout << "Rescaling old configuration to match desired temperature" << endl << endl;
        double sumv2 = 0.0, fs, xnew, ynew, znew;
        Move();
        for (int i=0; i<npart; ++i) sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
        sumv2 /= (double)npart;
        fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
        for (int i=0; i<npart; ++i){
          vx[i] *= fs;
          vy[i] *= fs;
          vz[i] *= fs;

          xnew = Pbc(x[i] - vx[i] * delta * 2);
          ynew = Pbc(y[i] - vy[i] * delta * 2);
          znew = Pbc(z[i] - vz[i] * delta * 2);

          x[i] = xold[i];
          y[i] = yold[i];
          z[i] = zold[i];

          xold[i] = xnew;
          yold[i] = ynew;
          zold[i] = znew;
        }
      }
    } else {
        cout << "No file old.0" << endl;
        old = false;
    }

  }

  //Prepare initial velocities
  if (old == false){
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rnd.Rannyu() - 0.5;
     vy[i] = rnd.Rannyu() - 0.5;
     vz[i] = rnd.Rannyu() - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }

  }

  rnd.SaveSeed();
  return;
}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew;
  std::vector<double> fx(npart), fy(npart), fz(npart);

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }

  return f;
}

void Measure(){ //Properties measurement
  double v, k, vij, wij, w;
  double dx, dy, dz, dr;

  v = 0.0; //reset observables
  k = 0.0;
  w = 0.0;
  //reset the hystogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     //update of the histogram of g(r)
     for(int bin=0; bin<nbins; ++bin) if ((dr >= bin*bin_size)&(dr<(bin+1)*bin_size)) walker[igofr+bin]+= 2;

     if(dr < rcut){
       vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
       wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);
      //Potential energy and virial
       v += vij;
       w += wij;
     }
    }
  }
//cycle over all particles
  for (int i=0; i<npart; ++i) k += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

//Measurement
  walker[iv] = 4.0*v; //Potential energy
  walker[ik] = k; //Kinetic energy
  walker[it] = (2.0 / 3.0) * k; //Temperature
  walker[ie] = walker[iv]+walker[ik]; //Total enery
  walker[iw] = 48.0 * w / 3.0; //Virial

  if (instant_print == true) {
    ofstream InstantTemp, InstantEtot;

    InstantEtot.open("output/instant.etot.0",ios::app);
    InstantTemp.open("output/instant.temp.0",ios::app);

    InstantEtot << walker[ie]/(double)npart << endl; //Total energy
    InstantTemp << walker[it]/(double)npart << endl; //Temp

    InstantEtot.close();
    InstantTemp.close();
  }
  return;
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
}

void Averages(int iblk) //Print results for current block
{
   ofstream Gofr, Gave, Epot, Pres, Temp, Etot, Ekin;
   const int wd=12;

    cout << "Block number " << iblk << endl;

    Epot.open("output/output.epot.0",ios::app);
    Pres.open("output/output.pres.0",ios::app);
    Gofr.open("output/output.gofr.0",ios::app);
    Gave.open("output/output.gave.0",ios::app);
    Etot.open("output/output.etot.0",ios::app);
    Ekin.open("output/output.ekin.0",ios::app);
    Temp.open("output/output.temp.0",ios::app);

    stima_epot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy
    glob_av[iv] += stima_epot;
    glob_av2[iv] += stima_epot*stima_epot;
    err[iv]=Error(glob_av[iv],glob_av2[iv],iblk);

    stima_ekin = blk_av[ik]/blk_norm/(double)npart; //Kinetic energy
    glob_av[ik] += stima_ekin;
    glob_av2[ik] += stima_ekin*stima_ekin;
    err[ik]=Error(glob_av[ik],glob_av2[ik],iblk);

    stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; //Pressure
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres*stima_pres;
    err[iw]=Error(glob_av[iw],glob_av2[iw],iblk);

    stima_etot = blk_av[ie]/blk_norm/(double)npart + vtail ;//Total energy
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot*stima_etot;
    err[ie]=Error(glob_av[ie],glob_av2[ie],iblk);

    stima_temp = blk_av[it]/blk_norm/(double)npart; //Temperature
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp*stima_temp;
    err[it]=Error(glob_av[it],glob_av2[it],iblk);

//Potential energy per particle
    Epot << setw(wd) << iblk <<  setw(wd) << stima_epot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err[iv] << endl;
//Pressure
    Pres << setw(wd) << iblk <<  setw(wd) << stima_pres << setw(wd) << glob_av[iw]/(double)iblk << setw(wd) << err[iw] << endl;
//Kinetic energy
    Ekin << setw(wd) << iblk <<  setw(wd) << stima_ekin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err[ik] << endl;
//Total energy
    Etot << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err[ie] << endl;
//Temperature
    Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err[it] << endl;

//g(r)
    double dv, r;
    r = bin_size;
    for (int bin=0; bin<nbins; ++bin){
      dv = 4*M_PI/3 *((r+bin_size)*(r+bin_size)*(r+bin_size)-r*r*r);
      stima_gofr = blk_av[igofr+bin]/blk_norm /(rho*npart*dv);
      Gofr << setw(wd) << bin <<  setw(wd) << r << setw(wd) << stima_gofr << endl;

      glob_av[igofr+bin] += stima_gofr;
      glob_av2[igofr+bin] += stima_gofr*stima_gofr;
      err[igofr+bin]=Error(glob_av[igofr+bin],glob_av2[igofr+bin],iblk);
      r += bin_size;
    }

//Final g(r)
    if (iblk == nblk){
      r = bin_size;
      for (int bin=0; bin<nbins; ++bin){
        Gave << setw(wd) << bin <<  setw(wd) << r << setw(wd) << glob_av[igofr+bin]/(double)iblk << setw(wd) << err[igofr+bin] << endl;
        r += bin_size;
      }
    }

    cout << "----------------------------" << endl << endl;

    Epot.close();
    Pres.close();
    Gofr.close();
    Gave.close();
    Etot.close();
    Ekin.close();
    Temp.close();
}

void Accumulate(void){ //Update block averages

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf,WriteOld;
  cout << endl;
  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  cout << "Print old configuration to file old.final " << endl << endl;
  WriteOld.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteOld <<  xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteOld.close();

  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("output/frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk){ //Uncertainty
 if (iblk==1) return 0;
 else return sqrt((sum2/(double)iblk - sum/(double)iblk*sum/(double)iblk)/(double)iblk);
}


/****************************************************************
*****************************************************************
*****************************************************************
*****************************************************************/
