/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019

Exercise 10.2
Parallel tempering for TSP
pt.cpp

Francesco Ariele Piziali
*****************************************************************
*****************************************************************/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "mpi.h"

#include "pt.hpp"

using namespace std;

/*****************************************************************/
int main(int argc, char* argv[]){

  MPI_Init(&argc,&argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &taskid );
  MPI_Comm_size( MPI_COMM_WORLD, &numtasks );

  Initialize();

  if ((nbeta != numtasks) and (taskid == 0)){//Exit if number of replicas not matching tasks
    cout << "ERROR: number of replicas to sample does not match the number of tasks requested" << endl;
    return 1;
  }

  for(int istep=1; istep <= nstep; ++istep){//Simulation
    MPI_Barrier(MPI_COMM_WORLD);

    tsp->Move();            //Metropolis move

    if(istep%nswap == 0) Swap(istep); //Swap attempt between replicas

    Averages(istep);        //Print results
  }


  Finalize();

  MPI_Finalize();
  return 0;
}
/*****************************************************************/

void Realization(){ //Realize an instance of the TSP
  Random rnd("seed.in", 0);
  if (geometry == 0){ //Cities on the unit circle
    double theta;
    for (int i = 0; i<ncities; ++i){
      theta = rnd.Rannyu(0,2*M_PI);
      cities_list[i] = City(sin(theta), cos(theta));
    }
  }
  else{ //Cities in a unit square centered in the origin
    for (int i = 0; i<ncities; ++i){
      cities_list[i] = City(rnd.Rannyu(-0.5,0.5), rnd.Rannyu(-0.5,0.5));
    }
  }
}

void Initialize(){//Read input, get a realization of TSP and create a population
  rnd = make_shared<Random>("seed.in", taskid+1);
  ifstream ReadInput;
  ReadInput.open("input.in");

  ReadInput >> ncities;
  ReadInput >> geometry;

  ReadInput >> nstep;

  ReadInput >> nbeta;

  beta.resize(nbeta);
  for(int i=0; i<nbeta; ++i) ReadInput >> beta[i];

  ReadInput >> nswap;

  ReadInput >> permutation;
  ReadInput >> inversion;
  ReadInput >> pshift;

  ReadInput.close();

  cities_list.resize(ncities);

  Realization();

  vector<double> rates({permutation, inversion, pshift});

  tsp = make_shared<Annealing>(rnd, cities_list, rates, beta[taskid]);

  Averages(0);

  accepted_swap = 0;
  attempted_swap = 0;
  min_length = tsp->Length();

  if(taskid == 0) Info();
  MPI_Barrier(MPI_COMM_WORLD);

  cout << "Process " << taskid + 1  << " of "  << numtasks << endl;
  cout << "Start - Length = " << tsp->Length() << endl << endl;
  MPI_Barrier(MPI_COMM_WORLD);
}

void Averages(int istep){//Print averages, min and max lengths to file; print configuration to file
  const int wd = 12;

  //Print to file
  ofstream Output;
  Output.open("output/output.len."+to_string(taskid+1),ios::app);
  Output  << setw(wd) << istep
          << setw(wd) << tsp->Beta()
          << setw(wd) << tsp->Length()
          << endl;
  Output.close();


  //Show during execution
  int iprint = nstep/100;
  if((istep%iprint == 0) and (istep!=0 ) and (taskid == 0)){
    cout  << "Step " << setw(6) << istep
          << " - Minimal length = " << setw(5) << min_length
          << endl;
    //tsp->Acceptance();
  }
}

void Finalize(){//Show results
  MPI_Barrier(MPI_COMM_WORLD);
  cout << "Process " << taskid + 1  << " of "  << numtasks << endl;
  cout << "Final - Length = " << tsp->Length() << endl << endl;

  MPI_Barrier(MPI_COMM_WORLD);
  if (taskid == 0){
    cout <<  "Minimal length found = " << setw(5) << min_length << endl;
    cout << "Swap acceptance rate = " << setw(5) << accepted_swap/attempted_swap << endl;
  }

}

void Info(){
  cout << "Travelling Salesman Problem - Parallel Tempering" << endl << endl;

  if(geometry==0) cout <<"Geometry: cities on a unit circle" << endl;
  else cout << "Geometry: cities inside a unit square" << endl;
  cout << "Number of cities to visit = " << ncities << endl << endl;

  cout << "Number of Metropolis steps  = " << nstep << endl;
  cout << "Number of replicas = " << nbeta << endl;
  cout << "Replicas:" << endl;
  for (int i=0; i<nbeta; ++i) cout << "   " << i + 1 << ": beta = " << beta[i] << endl;
  cout << "Swap attempt between replicas every " << nswap << " steps" << endl;
  cout << endl;

  cout << "Rates:" << endl;
  cout << "pair permutation = " << permutation << endl;
  cout << "inversion = " <<  inversion << endl;
  cout << "partial shift = " <<  pshift << endl;
  cout << endl;
}

void Swap(int istep){
  int root = 0;

  std::vector<int> task_config = tsp->GetConfiguration();
  double task_length = tsp->Length();

  std::vector<int>  gather_config(ncities*numtasks);
  std::vector<double> gather_lenght(numtasks);

  MPI_Gather( &task_length, 1, MPI_DOUBLE, gather_lenght.data(), 1, MPI_DOUBLE, root, MPI_COMM_WORLD );
  MPI_Gather( &task_config[0], ncities, MPI_INTEGER, &gather_config[0], ncities, MPI_INTEGER, root, MPI_COMM_WORLD );

  if(taskid == root){
    int r1 = rnd->Integer(1, numtasks-1); //recieve configuration (<T, >beta)
    int r2 = r1 - 1; //send configuration (>T, <beta)
    double d_beta = beta[r1]-beta[r2];
    double d_energy = gather_lenght[r1]-gather_lenght[r2];
    double p = exp(d_beta*d_energy);

    //for(int i=0; i<gather_lenght.size(); ++i) cout << gather_lenght[i] << "  ";
    //cout << endl;

    double temp_min = *min_element(gather_lenght.begin(), gather_lenght.end());

    if (temp_min< min_length){//save minimal lenght and best config
      ofstream MinOut, Best;

      min_length = temp_min;
      int min_pos = distance(gather_lenght.begin(), min_element(gather_lenght.begin(), gather_lenght.end()));

      MinOut.open("output/output.len.0", ios::app);
      MinOut << setw(12) << istep << setw(6) << min_pos+1 <<setw(12) << min_length <<  endl;
      MinOut.close();

      std::vector<int>  best_config(ncities);
      copy(gather_config.begin()+min_pos*ncities, gather_config.begin()+min_pos*ncities +ncities, best_config.begin());
      Route best(cities_list, best_config);
      best.Print(min_pos+1, istep);
    }

    if (rnd->Rannyu() <= p) {
      //cout << r1 << "  " << r2 << "   "<< p << endl;
      std::vector<int> temp1(ncities), temp2(ncities);
      for(int i=0; i<ncities; ++i){//copy config
        temp1[i] = gather_config[r2*ncities + i];
        temp2[i] = gather_config[r1*ncities + i];
      }
      for(int i=0 ; i<ncities; ++i) {//paste config
        gather_config[r1*ncities + i] = temp1[i];
        gather_config[r2*ncities + i] = temp2[i];
      }

      ofstream SwapOut;
      SwapOut.open("output/swap.log",ios::app);
      //cout << "Configuration swapped from " << r2 << " to " << r1 << "!" << endl;
      SwapOut << setw(12) << istep << setw(6) << r2+1 <<setw(6) << r1+1 <<  endl;
      SwapOut.close();
      accepted_swap ++;
    }
    attempted_swap++;
  }

  MPI_Scatter(&gather_config[0], ncities, MPI_INTEGER, &task_config[0], ncities, MPI_INTEGER, root, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  tsp->SetConfiguration(task_config);
  MPI_Barrier(MPI_COMM_WORLD);
}
