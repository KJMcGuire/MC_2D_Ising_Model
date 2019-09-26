///Author:  Kellie McGuire     kellie@kelliejensen.com

///This program computes the energy, magnetization, magnetic susceptibility, and specific heat
///for a 2D Ising Model of a ferromagnet using Monte Carlo methods with
///Metropolis sampling. Code was adapted from "Introduction to Monte Carlo
///methods for an Ising Model of a Ferromagnet'' by Jacques Kotze.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>   //Need for rand; srand
#include <iostream>
#include <fstream>    //Need for ofstream
#include <time.h>     //Need for time





///Structure for 2D lattice
struct lattice_type {
  int x;
  int y;
};


///Define simulation parameters
const int size = 128;    //Lattice size;
int n = size*size;     //total number of lattice points
int lat[size+1][size+1];
float maxT=5.0;        //Max temp
float minT=0.5;        //Min temp
float Tchange=0.1;     //Step size for temp
long int mcs=10000;    //Number of Monte Carlo steps
int skip=1000;         //Number of steps to omit from data (thermalizing)
float j=-1;             //Coupling constant  (>0 for ferromagnetic; <0 for antiferro)
float de = 0.0;        //Initialize change in energy
double norm = 1.0/float(mcs*n);


double MyRand(){
  double r = rand()/(RAND_MAX+1.0);
  return r;
}

///Initialize lattice with random spins and print to screen
void init_hot(int lat[size+1][size+1]){
  int i, j;
    for (j=size; j>=1; j--){    //Column-major loop
      for (i=1; i<=size; i++){
        if(MyRand() < 0.5)
            lat[i][j]=1;
            else
              lat[i][j]=-1;
        printf("%2d ", lat[i][j]);
            }
        printf("\n");
    }
}

///Initialize lattize w/ all spins at 1
void init_cold(int lat[size+1][size+1]){
  int i, j;
    for (j=size; j>=1; j--){    //Column-major loop
      for (i=1; i<=size; i++){
            lat[i][j]=1;
        printf("%2d ", lat[i][j]);
            }
        printf("\n");
    }
}

///Pick random position
void random_position(lattice_type &pos){
  pos.x=(int)ceil(MyRand()*size);
  pos.y=(int)ceil(MyRand()*size);
  if(pos.x>size || pos.y>size){
    printf("Point falls outside array.");
    exit;
  }
}

///Calculate energy of lattice position
int energy_pos(lattice_type &pos){
  int left, right, lower, upper;
  ///Periodic boundary conditions
  if(pos.y==size) upper=1;
    else upper=pos.y+1;
  if(pos.y==1) lower=size;
    else lower=pos.y-1;
  if(pos.x==size) right=1;
    else right=pos.x+1;
  if(pos.x==1) left=size;
    else left=pos.x-1;

  ///six nearest neighbors for triangle lattice
float e = -j*lat[pos.x][pos.y]*(lat[pos.x][upper]+//
          lat[right][upper]+lat[left][lower]+//
          lat[pos.x][lower]+lat[right][pos.y]+lat[left][pos.y]);



  ///Four nearest neighbors for square lattice
// float e = -j*lat[pos.x][pos.y]*(lat[pos.x][upper]+lat[pos.x][lower]//
//   +lat[right][pos.y]+lat[left][pos.y]);


 return e;
}


///Calculate probability of flipping
bool p_flip(lattice_type pos, float &de, float &T){
  de = -2*j*energy_pos(pos);   //Change in energy for spin
  if(de < 0) return true;    //Flip if delta E less than 0
    else if(MyRand()<exp(-de/T)) return true;   //Flip due to heat bath
      else return false;
}


///Flip lattice value
void flip(lattice_type pos){
  lat[pos.x][pos.y]=-lat[pos.x][pos.y];
}


///Calculate magnetization
int magnetization(){
  int m = 0;
  for(int y=size; y>= 1; y--){
    for(int x=1; x<= size; x++){
      m += lat[x][y];
    }
  }
  return m;
}

///Calculate square of magnetization
int sq_mag(){
  int m = 0;
  for(int y=size; y>= 1; y--){
    for(int x=1; x<= size; x++){
      m += lat[x][y];
    }
  }
  return m*m;
}


///Calculate lattice energy
int total_energy(){
  lattice_type pos;
  int E=0;
  for(int y = size; y >= 1; y--){
    pos.y=y;
    for(int x=1; x <= size; x++){
      pos.x=x;
      E+=energy_pos(pos);
    }
  }
return E;
}


///Main function
int main(int argc, char *argv[]){

if (argc < 2){
  printf("Specify an output file name\n");
  exit;
}

///Create output file for data
std::ofstream outfile (argv[1]);

///Create variables for storing changes in observables
double Mag=0, Mag_avg=0, E_avg_sq=0, EsqAvg=0, C=0, E=0;
double MsqAvg=0, Chi=0, E_avg=0, M_avg_sq=0, Mag_abs=0;

///Initialize lattice
srand(time(NULL));  //Rand() seed
lattice_type pos;
init_hot(lat);
//init_cold(lat);




///Temp loop
for(float T=maxT; T>=minT; T=T-Tchange){

  ///Allow lattice to thermalize
  for(int a=1; a<=skip; a++){
    for(int b=1; b<n; b++){
      random_position(pos);
      if(p_flip(pos, de, T)) flip(pos);   //Flip spin at lattice point
    }
  }

    ///Initialize variables for storing changes in observables
    int Mag=0;
    double AvgMag=0, Energy=0, Sq_Energy=0, Sq_Mag=0, Mag_abs=0;
    E= total_energy();

  ///Monte Carlo loop
  for(int i=1; i<=mcs; i++){

    ///Metropolis loop
    for(int i=1; i<=n; i++)
    {
      random_position(pos);
      if(p_flip(pos, de, T))
      {
       flip(pos);   //Flip spin at lattice point
       E += 2*de;
      }
    }
    Mag += magnetization();
    Mag_abs += abs(magnetization());
    Energy += E/2.0; //To avoid double counting
    Sq_Energy += E*E/4.0;
    Sq_Mag += sq_mag();


  }

  Mag_avg = Mag_abs*norm;
  E_avg = Energy*norm;
  E_avg_sq = E_avg*E_avg;
  M_avg_sq = Mag_avg*Mag_avg;
  EsqAvg = Sq_Energy*norm;
  MsqAvg = Sq_Mag*norm*norm*mcs;
  C = (EsqAvg-(E_avg_sq*n))/(T*T);  //Specific heat per spin
  printf("%2f ", C);
  Chi = ((MsqAvg-(Mag_avg*Mag_avg))*n)/T; //Magnetic susceptibility

  ///Write data to file
  outfile<<T<<"\t"<<Mag_avg<<"\t"<<E_avg<<"\t"<<C<<"\t"<<
  Chi<<"\t"<<std::endl;
}
outfile.close();
return 0;
}
