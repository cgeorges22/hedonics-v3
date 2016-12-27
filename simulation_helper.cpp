
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>  //JKR 10/7/16
#include <sstream> //JKR 10/7/16
#include "hedonics-firm-v2.10.h"
#include "hedonics-utilities-v2.10.h"
#include "time.h"
#include "mpi.h"
#include <cassert> //7/4/09
#include <cstring>
#include <deque> //JKR 10/31/16
#include <cmath>
#include "master_slave_hedonics.h"
#define MASTER 0
#define WORKTAG 2
#define DIETAG 3
#define NEEDWORK 4
using namespace std;


void master(int arg1, int arg2);
void slave(int rank);
void getSeeds(int& start, int& end);

int main(int argc, char *argv[]){
  
  int size, rank, perJobs;
  int start, end;

  getSeeds(start,end);
  //Got the inputs

  MPI_Init(0, 0);  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
 
  if(rank == MASTER) master(start, end);
  else slave(rank);
 
  MPI_Finalize();
  return 0;

}

void master(int start, int end){
  int totalJobs, rank, ntasks;
  deque<int>globalQueue;
  totalJobs = abs(end-start);

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
 

  // Fills the Global Queue
  int temp = start;                                    
  for(int i = 0; i <= totalJobs; i++, temp++){
    globalQueue.push_back(temp); 
  } 
    
  // Intialize mass seed distribution
  for(rank = 1; rank < ntasks; rank++){
    if(!globalQueue.empty()){
      int buffer = globalQueue.front();
      globalQueue.pop_front();
      MPI_Send(&buffer, 1, MPI_INT,rank, 1, MPI_COMM_WORLD);
    }
  }

  // If there are more seeds than processors, continue to distribute
  while(!globalQueue.empty()){ 
    int buffer = globalQueue.front();
    globalQueue.pop_front();
    MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, NEEDWORK,MPI_COMM_WORLD,&status);

    MPI_Send(&buffer, 1, MPI_INT, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD); 
  }

  // Kill all other processors once all simulations are over
  for(rank = 1; rank < ntasks; rank++ ){
    MPI_Send(0,0,MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
  } 

  
    
   return;
  } //end of master function


void slave(int rank){

   int randSeed;
   MPI_Status status;

   for(;;){

     // Receive RandSeed for simulation
     MPI_Recv(&randSeed,1,MPI_INT,MASTER,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      
     if(status.MPI_TAG == DIETAG){
       return;
     }
     simulation(randSeed, rank);

     // Send back to Master that simulation is over
     MPI_Send(0,0,MPI_INT,MASTER,NEEDWORK,MPI_COMM_WORLD);
   }
}


void getSeeds(int& start, int& end){
  char line[256];
  char* variable;
  char* value;
  // Uses the file "input.txt"
  std::fstream input;
  input.open("input.txt", std::fstream::in);
  
  // Go through the input file line by line, tokenize the line into three strings
  // Taken from hedonics-v2.10
  input.getline(line, 256);
  variable = strtok(line, " ");

  while (strcmp(variable, "randSeedStart")) {
    input.getline(line, 256);
    variable = strtok(line, " ");
  }
  start = atoi(strtok(NULL, " "));
  input.getline(line, 256);
  strtok(line, " ");
  end = atoi(strtok(NULL, " ")); 
}

