/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019
Excercise 04.1
Molecular_Dynamics_NVE.h
Francesco Ariele Piziali
*****************************************************************
*****************************************************************/

#include <vector>
//parameters, observables
//const int m_props=1000;
int n_props;
int iv,ik,it,ie, iw, igofr, nbins;
double stima_epot, stima_ekin, stima_etot, stima_temp, stima_pres, stima_gofr;
double vtail,ptail,bin_size;
std::vector<double>  walker;

// averages
std::vector<double> blk_av;
double blk_norm;
std::vector<double> glob_av,glob_av2, err;

//configuration
//const int m_part=108;
std::vector<double> x, y, z, xold, yold, zold;
std::vector<double> vx, vy, vz;

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, nblk, iprint;
double delta;
bool old,rescale,print_xyz, instant_print;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);

//statistics
void Reset(int);
void Accumulate(void);
void Averages(int);
double Error(double,double,int);

/****************************************************************
*****************************************************************
*****************************************************************
*****************************************************************/
