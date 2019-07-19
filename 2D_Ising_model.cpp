//Author:  Kellie McGuire     kellie@kelliejensen.com

//This program computes the energy, magnetic susceptibility, and specific heat
//for a 2D Ising Model of a ferromagnet using Monte Carlo methods with
//Metropolis sampling. Code was adapted from "Introduction to Monte Carlo
//methods for an Ising Model of a Ferromagnet'' by Jacques Kotze.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>   //Need for rand; srand
#include <iostream>
#include <fstream>    //Need for ofstream
#include <time.h>     //Need for time


//Structure for 2D lattice
struct lattice_type {
  int x;
  int y;
};


//Define simulation parameters
const int size = 32;    //Lattice size;
int n = size*size;     //total number of lattice points
int lat[size+1][size+1];
float maxT=5.0;        //Max temp
float minT=0.5;        //Min temp
float Tchange=0.1;     //Step size for temp
long int mcs=10000;    //Number of Monte Carlo steps
int skip=1000;         //Number of steps to omit from data (thermalize)
float j=1;             //Coupling constant
float de = 0.0;        //Initialize change in energy
double norm = 1.0/float(mcs*n);

double MyRand(){
  double r = rand()/(RAND_MAX+1.0);
  return r;
}

//Initialize lattice with random input and print to screen
void initialize(int lat[size+1][size+1]){
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

//Establish random position
void random_position(lattice_type &pos){
  pos.x=(int)ceil(MyRand()*size);
  pos.y=(int)ceil(MyRand()*size);
  if(pos.x>size || pos.y>size){
    printf("Starting position falls outside array.");
    exit;
  }
//printf("%d %d", pos.x,pos.y);
}

//Calculate energy of lattice position
int energy_pos(lattice_type &pos){
  int left, right, lower, upper;
  //Periodic boundary conditions
  if(pos.y==size) upper=1;
    else upper=pos.y+1;
  if(pos.y==1) lower=size;
    else lower=pos.y-1;
  if(pos.x==size) right=1;
    else right=pos.x+1;
  if(pos.x==1) left=size;
    else left=pos.x-1;

  float e = -j*lat[pos.x][pos.y]*(lat[pos.x][upper]+lat[pos.x][lower]//
    +lat[right][pos.y]+lat[left][pos.y]);
  return e;
}

//Calculate probability of flipping
bool p_flip(lattice_type pos, float &de, float &T){
  de = -2*energy_pos(pos);   //Change in energy for spin
  if(de < 0) return true;    //Flip if delta E less than 0
    else if(MyRand()<exp(-de/T)) return true;   //Flip due to heat bath
      else return false;
}


//Flip lattice value
void flip(lattice_type pos){
  lat[pos.x][pos.y]=-lat[pos.x][pos.y];
}


//Calculate magnetization
int magnetization(){
  int m = 0;
  for(int y=size; y>= 1; y--){
    for(int x=1; x<= size; x++){
      m += lat[x][y];
    }
  }
  return m;
}



//Main function
int main(int argc, char *argv[]){

if (argc < 2){
  printf("Specify an output file name\n");
  exit;
}

//Create output file for data
std::ofstream outfile (argv[1]);

//Initialize lattice
srand(time(NULL));  //Random seed
lattice_type pos;
initialize(lat);

//Allow lattice to thermalize
for(int i=1; i<skip; i++){
  random_position(pos);
  if(p_flip(pos, de, maxT)) flip(pos);   //Flip spin at lattice point
  }

//Temp loop
for(float T=maxT; T>=minT; T=T-Tchange){

  //Create variables for storing changes in observables
  int Mag=0;
  float AvgMag=0.0;

  //Monte Carlo loop
  for(int i=1; i<mcs; i++){

    //Metropolis loop
    for(int i=1; i<n; i++){
      random_position(pos);
      if(p_flip(pos, de, T)) flip(pos);   //Flip spin at lattice point
      }
    Mag += magnetization();
  }
  AvgMag = Mag*norm;


  //Write data to file
  outfile << T <<"\t"<<abs(AvgMag)<< std::endl;
}
outfile.close();
return 0;
}
