/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019
Excercise 05.1
main.cpp
Francesco Ariele Piziali
*****************************************************************
*****************************************************************/

#include "main.h" //Global variables definition


int main (int argc, char *argv[]){

//Input parameters
  if (argc < 4) {
    cerr << "USAGE: " << argv[0] << " [#step] [delta 1s] [delta 2p] [delta 3d]"<< endl;
    return 1;
  }

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

  M = atoi(argv[1]);
  r.resize(M); //Number of Metropolis steps
  AccOut.open("acceptance.out");

  cout << "Metropolis estimation of <r> for 1s and 2p Hydrogen atom orbitals" << endl;
  cout << "Number of Metropolis steps = " << M << endl << endl;

  AccOut << M << endl;

  //1s
  delta = atof(argv[2]); //Step lenght;

  CoorOut.open("coordinates_1s.xyz");
  cout << "Uniform transition probability step = " << delta << endl;
  AccOut << delta << endl;

  accepted = 0;
  x_old = 0;
  y_old = 0;
  z_old = 1;

  p_old = p_1s(x_old, y_old, z_old);

  r[0] = distance(x_old, y_old, z_old);

  for (int i=1; i<M; ++i){
    //Random solid angle
    u = rnd.Rannyu(-1,1);
    theta = rnd.Rannyu(0,2*M_PI);

    x_new = x_old + delta * sqrt(1-u*u) * cos(theta);
    y_new = y_old + delta * sqrt(1-u*u) * sin(theta);
    z_new = z_old + delta * u;

    p_new = p_1s(x_new, y_new, z_new);

    accept = test(p_new, p_old, rnd.Rannyu());

    if (accept == true){
      x_old = x_new;
      y_old = y_new;
      z_old = z_new;

      r[i] = distance(x_new, y_new, z_new);
      p_old = p_new;

      accepted++;
    } else {
      r[i] = r[i-1];
    }
    if(i%100 == 0)  CoorOut << x_old << "\t" << y_old << "\t" << z_old << endl;
  }

  CoorOut.close();
  cout << "1s acceptance rate = " << accepted/(double)M << endl;
  AccOut << accepted/(double)M << endl;

  blocking(r, N, "r_1s.out");

  cout << endl;

  //2p
  delta = atof(argv[3]); //Step lenght;

  CoorOut.open("coordinates_2p.xyz");
  cout << "Uniform transition probability step = " << delta << endl;
  AccOut << delta << endl;

  accepted = 0;

  x_old = 0;
  y_old = 0;
  z_old = 3;

  p_old = p_2p(x_old, y_old, z_old);

  r[0] = distance(x_old, y_old, z_old);

  for (int i=1; i<M; ++i){
    //Random solid angle
    u = rnd.Rannyu(-1,1);
    theta = rnd.Rannyu(0,2*M_PI);

    x_new = x_old + delta * sqrt(1-u*u) * cos(theta);
    y_new = y_old + delta * sqrt(1-u*u) * sin(theta);
    z_new = z_old + delta * u;

    p_new = p_2p(x_new, y_new, z_new);

    accept = test(p_new, p_old, rnd.Rannyu());

    if (accept == true){
      x_old = x_new;
      y_old = y_new;
      z_old = z_new;

      r[i] = distance(x_new, y_new, z_new);
      p_old = p_new;

      accepted++;
    } else {
      r[i] = r[i-1];
    }
    if(i%100 == 0)  CoorOut << x_old << "\t" << y_old << "\t" << z_old << endl;
  }

  CoorOut.close();
  cout << "2p acceptance rate = " << accepted/(double)M << endl;
  AccOut << accepted/(double)M << endl;

  blocking(r, N, "r_2p.out");

  cout << endl;

  //1s normal transition
  delta = atof(argv[2])/2; //Step lenght;
  double r_ave = 0;

  CoorOut.open("coordinates_1s_norm.xyz");
  cout << "Normal transition probability N(0," << delta << ")" << endl;
  AccOut << delta << endl;

  accepted = 0;
  x_old = 0;
  y_old = 0;
  z_old = 1;

  p_old = p_1s(x_old, y_old, z_old);

  r[0] = distance(x_old, y_old, z_old);

  for (int i=1; i<M; ++i){
    //Normal transition for each coordinate
    tx = rnd.Gauss(0,delta);
    ty = rnd.Gauss(0,delta);
    tz = rnd.Gauss(0,delta);
    x_new = x_old + tx;
    y_new = y_old + ty;
    z_new = z_old + tz;

    p_new = p_1s(x_new, y_new, z_new);

    accept = test(p_new, p_old, rnd.Rannyu());
    r_ave += distance(x_new-x_old, y_new-y_old, z_new-z_old);
    if (accept == true){

      x_old = x_new;
      y_old = y_new;
      z_old = z_new;

      r[i] = distance(x_new, y_new, z_new);
      p_old = p_new;

      accepted++;
    } else {
      r[i] = r[i-1];
    }
    if(i%100 == 0)  CoorOut << x_old << "\t" << y_old << "\t" << z_old << endl;
  }

  CoorOut.close();
  cout << "1s acceptance rate  = " << accepted/(double)M << endl;
  cout << "Average step lenght = " << r_ave/(double)M << endl;
  AccOut << accepted/(double)M << endl;
  AccOut << r_ave/(double)M << endl;
  blocking(r, N, "r_1s_norm.out");

  cout << endl;

  //2p normal transition
  delta = atof(argv[3])/2; //Step lenght;

  CoorOut.open("coordinates_2p_norm.xyz");
  cout << "Normal transition probability N(0," << delta << ")" << endl;
  AccOut << delta << endl;

  accepted = 0;
  x_old = 0;
  y_old = 0;
  z_old = 3;
  r_ave = 0;
  p_old = p_2p(x_old, y_old, z_old);

  r[0] = distance(x_old, y_old, z_old);

  for (int i=1; i<M; ++i){
    //Normal transition for each coordinate
    tx = rnd.Gauss(0,delta);
    ty = rnd.Gauss(0,delta);
    tz = rnd.Gauss(0,delta);
    x_new = x_old + tx;
    y_new = y_old + ty;
    z_new = z_old + tz;

    p_new = p_2p(x_new, y_new, z_new);

    accept = test(p_new, p_old, rnd.Rannyu());
    r_ave += distance(x_new-x_old, y_new-y_old, z_new-z_old);
    if (accept == true){
      x_old = x_new;
      y_old = y_new;
      z_old = z_new;

      r[i] = distance(x_new, y_new, z_new);
      p_old = p_new;

      accepted++;
    } else {
      r[i] = r[i-1];
    }
    if(i%100 == 0)  CoorOut << x_old << "\t" << y_old << "\t" << z_old << endl;
  }

  CoorOut.close();
  cout << "2p acceptance rate = " << accepted/(double)M << endl;
  cout << "Average step lenght = " << r_ave/(double)M << endl;
  AccOut << accepted/(double)M << endl;
  AccOut << r_ave/(double)M << endl;
  blocking(r, N, "r_2p_norm.out");

  cout << endl;

  //3d 3,2,0 uniform
  delta = atof(argv[4]); //Step lenght;

  CoorOut.open("coordinates_320.xyz");
  cout << "Uniform transition probability step = " << delta << endl;
  AccOut << delta << endl;

  accepted = 0;

  x_old = 0;
  y_old = 0;
  z_old = 3;

  p_old = p_320(x_old, y_old, z_old);

  r[0] = distance(x_old, y_old, z_old);

  for (int i=1; i<M; ++i){
    //Random solid angle
    u = rnd.Rannyu(-1,1);
    theta = rnd.Rannyu(0,2*M_PI);

    x_new = x_old + delta * sqrt(1-u*u) * cos(theta);
    y_new = y_old + delta * sqrt(1-u*u) * sin(theta);
    z_new = z_old + delta * u;

    p_new = p_320(x_new, y_new, z_new);

    accept = test(p_new, p_old, rnd.Rannyu());

    if (accept == true){
      x_old = x_new;
      y_old = y_new;
      z_old = z_new;

      r[i] = distance(x_new, y_new, z_new);
      p_old = p_new;

      accepted++;
    } else {
      r[i] = r[i-1];
    }
    if(i%100 == 0)  CoorOut << x_old << "\t" << y_old << "\t" << z_old << endl;
  }

  CoorOut.close();
  cout << "3d (3,2,0) acceptance rate = " << accepted/(double)M << endl;
  AccOut << accepted/(double)M << endl;

  blocking(r, N, "r_320.out");
  cout << endl;

  //3d 3,2,1 uniform
  delta = atof(argv[4]); //Step lenght;

  CoorOut.open("coordinates_321.xyz");
  cout << "Uniform transition probability step = " << delta << endl;
  AccOut << delta << endl;

  accepted = 0;

  x_old = 0;
  y_old = 0;
  z_old = 3;

  p_old = p_321(x_old, y_old, z_old);

  r[0] = distance(x_old, y_old, z_old);

  for (int i=1; i<M; ++i){
    //Random solid angle
    u = rnd.Rannyu(-1,1);
    theta = rnd.Rannyu(0,2*M_PI);

    x_new = x_old + delta * sqrt(1-u*u) * cos(theta);
    y_new = y_old + delta * sqrt(1-u*u) * sin(theta);
    z_new = z_old + delta * u;

    p_new = p_321(x_new, y_new, z_new);

    accept = test(p_new, p_old, rnd.Rannyu());

    if (accept == true){
      x_old = x_new;
      y_old = y_new;
      z_old = z_new;

      r[i] = distance(x_new, y_new, z_new);
      p_old = p_new;

      accepted++;
    } else {
      r[i] = r[i-1];
    }
    if(i%100 == 0)  CoorOut << x_old << "\t" << y_old << "\t" << z_old << endl;
  }

  CoorOut.close();
  cout << "3d (3,2,1) acceptance rate = " << accepted/(double)M << endl;
  AccOut << accepted/(double)M << endl;

  blocking(r, N, "r_321.out");

  rnd.SaveSeed();
  return 0;
}

/****************************************************************/
//FUNCTIONS

double error(double sum, double sum2, int n){
 if (n==0) return 0;
 else return (sqrt((sum2 - sum*sum)/n));
}

double distance(double x, double y, double z){
  return sqrt(x*x + y*y + z*z);
}

bool test(double p_new, double p_old, double t){
  bool result = false;
  if (t<fmin(1, p_new/p_old)) result = true;
  return result;
}

double p_1s(double x, double y, double z){
  double r = distance(x,y,z);
  return exp(-2*r)/M_PI;
}

double p_2p(double x, double y, double z){
  double r = distance(x,y,z);
  double theta = atan(y/x);
  return (r*r*exp(-r)*cos(theta)*cos(theta)/32/M_PI);
}

double p_320(double x, double y, double z){
  double r = distance(x,y,z);
  double r2 = r*r;
  double theta = atan(y/x);
  return r2*r2*exp(-2*r/3)*(3*cos(theta)*cos(theta)-1)*(3*cos(theta)*cos(theta)-1)/(81*81*6*M_PI);
  //return r2*r2*exp(-2*r/3)*sin(theta)*sin(theta)*cos(theta)*cos(theta)/(81*81*M_PI);
}
double p_321(double x, double y, double z){
  double r = distance(x,y,z);
  double r2 = r*r;
  double theta = atan(y/x);
  //return r2*r2*exp(-2*r/3)*(3*cos(theta)*cos(theta)-1)*(3*cos(theta)*cos(theta)-1)/(81*81*6*M_PI);
  return r2*r2*exp(-2*r/3)*sin(theta)*sin(theta)*cos(theta)*cos(theta)/(81*81*M_PI);
}

void blocking(std::vector<double> var, int nblocks, std::string output_file){
  std::vector<double> sum_prog(nblocks), sum2_prog(nblocks), err_prog(nblocks), ave_blocks(nblocks);
  double sum;
  int nelements = var.size()/nblocks;
  int k;

  ofstream output;
  output.open(output_file);

  for (int i=0; i<nblocks; ++i){
    sum = 0;
    for (int j=0; j<nelements; ++j){
      k = j + i*nelements;
      sum += var[k];
    }
    ave_blocks[i] = sum/double(nelements);
   }

  for (int i=0; i<nblocks; ++i){
    for (int j=0; j<(i+1); ++j){
      sum_prog[i] += ave_blocks[j];
      sum2_prog[i] += ave_blocks[j]*ave_blocks[j];
    }
      sum_prog[i] /= (i+1);
      sum2_prog[i] /= (i+1);
      err_prog[i] = error(sum_prog[i], sum2_prog[i], i);
      output << sum_prog[i] << "\t" << err_prog[i] << endl;
  }
  output.close();

  cout << "Blocking statistics saved in " << output_file << endl;

  return;
}



/****************************************************************
*****************************************************************
*****************************************************************
*****************************************************************/
