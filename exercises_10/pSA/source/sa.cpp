/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019

Exercise 10.1
Simulated annealing for TSP
sa.cpp

Francesco Ariele Piziali
*****************************************************************
*****************************************************************/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "mpi.h"

#include "sa.hpp"

using namespace std;

/*****************************************************************/
int main(int argc, char* argv[]){

  MPI_Init(&argc,&argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &taskid );
  MPI_Comm_size( MPI_COMM_WORLD, &numtasks );

  Initialize();

  for(int istep=1; istep <= nstep; ++istep){//Simulation

    if (istep%ncool == 0){//Cooling
      MPI_Barrier(MPI_COMM_WORLD);
      beta = beta / cool;    //pow(0.5, istep/ncool);
      tsp->SetBeta(beta);
      Gather(istep);
    }

    tsp->Move();            //Metropolis move

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
  ReadInput >> beta;
  ReadInput >> cool;
  ReadInput >> ncool;

  ReadInput >> permutation;
  ReadInput >> inversion;
  ReadInput >> pshift;

  ReadInput.close();

  cities_list.resize(ncities);

  Realization();

  vector<double> rates({permutation, inversion, pshift});

  tsp = make_shared<Annealing>(rnd, cities_list, rates, beta);

  Averages(0);
  min_length = tsp->Length();


  if(taskid == 0) Info();
  MPI::COMM_WORLD.Barrier();

  cout << "Process " << taskid + 1  << " of "  << numtasks << endl;
  cout << "Start - Length = " << tsp->Length() << endl << endl;
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

  int iprint = nstep/10;
  if((istep%iprint == 0) and (istep!=0 ) and (taskid == 0)){
    cout  << "Step " << setw(6) << istep
          << " - Minimal length = " << setw(5) << min_length
          << endl;
  }

  if(istep == nstep) tsp->Print(taskid+1,istep);
}

void Info(){
  cout << "Travelling Salesman Problem - Parallel Simulated Annealing" << endl << endl;

  if(geometry==0) cout <<"Geometry: cities on a unit circle" << endl;
  else cout << "Geometry: cities inside a unit square" << endl;
  cout << "Number of cities to visit = " << ncities << endl;

  cout << "Number of Metropolis steps  = " << nstep << endl;
  cout << "Initial beta = " << beta << endl;
  cout << "Cooling rate = " << cool << " every " << ncool << " steps " << endl << endl;

  cout << "Rates:" << endl;
  cout << "pair permutation = " << permutation << endl;
  cout << "inversion = " <<  inversion << endl;
  cout << "partial shift = " <<  pshift << endl;
  cout << endl;
}

void Finalize(){//Show results
  MPI_Barrier(MPI_COMM_WORLD);
  cout << "Process " << taskid + 1  << " of "  << numtasks << endl;
  cout << "Final - Length = " << tsp->Length() << endl << endl;

  MPI_Barrier(MPI_COMM_WORLD);
  if (taskid == 0){
    cout << endl;
    cout <<  "Minimal length found = " << setw(5) << min_length << endl;
  }
}

void Gather(int istep){
  int root = 0;

  double task_length = tsp->Length();
  std::vector<int> task_config = tsp->GetConfiguration();

  std::vector<double> gather_lenght(numtasks);
  std::vector<int>  gather_config(ncities*numtasks);

  MPI_Gather( &task_length, 1, MPI_DOUBLE, gather_lenght.data(), 1, MPI_DOUBLE, root, MPI_COMM_WORLD );
  MPI_Gather( &task_config[0], ncities, MPI_INTEGER, &gather_config[0], ncities, MPI_INTEGER, root, MPI_COMM_WORLD );

  if(taskid == root){
    double temp_min = *min_element(gather_lenght.begin(), gather_lenght.end());

    if (temp_min< min_length){
      ofstream MinOut;
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
  }

}
