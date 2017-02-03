
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


void master(int arg1, int arg2, double paramStart, double paramStop, int paramNum);
void slave(int rank);
void getSeeds(int& start, int& end, string& paramName, double& paramStart, double& paramStop, int& paramNum);
//MPI_File file;
MPI_File * files; 
void openFiles(double paramStart, double paramStop, int paramNum);
void closeFiles(int paramNum);

int main(int argc, char *argv[]){
  
  int size, rank, perJobs;
  int start, end;
  string paramName;
  double paramStart, paramStop;
  int paramNum;

  getSeeds(start, end, paramName, paramStart, paramStop, paramNum);
  //Got the inputs

  MPI_Init(0, 0);  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // open the files - must be done by all procesoors synchronously
  files = new MPI_File [paramNum];
  openFiles(paramStart, paramStop, paramNum);

  if(rank == MASTER) master(start, end, paramStart, paramStop, paramNum);
  else slave(rank);

  closeFiles(paramNum); 

  MPI_Finalize();
  return 0;

}

void master(int start, int end, double paramStart, double paramStop, int paramNum){
  int totalJobs, rank, ntasks, paramSpace;
  deque<int>seedQueue;
  deque<double>paramQueue;
  totalJobs = abs(end-start);

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
 

  // Fills the Global Queue
  int seed = start;                                   
  double param; 
  double step = (paramStop - paramStart)/(paramNum - 1);

  for(int i = 0; i <= totalJobs; i++, seed++){
    param = paramStart;
    for (int j = 0; j < paramNum; j++, param += step) {
      seedQueue.push_back(seed); 
      seedQueue.push_back(j); 
      paramQueue.push_back(param);
    }
  } 

  int fileNum;
  // Intialize mass seed distribution
  for(rank = 1; rank < ntasks; rank++){
    if(!seedQueue.empty()){
      seed = seedQueue.front();
      seedQueue.pop_front();

      fileNum = seedQueue.front();
      seedQueue.pop_front();

      param = paramQueue.front();
      paramQueue.pop_front();
   
      MPI_Send(&seed, 1, MPI_INT,rank, 1, MPI_COMM_WORLD);
      MPI_Send(&fileNum, 1, MPI_INT,rank, 1, MPI_COMM_WORLD);
      MPI_Send(&param, 1, MPI_DOUBLE,rank, 1, MPI_COMM_WORLD); 
    }
  }

  // If there are more seeds than processors, continue to distribute
  while(!seedQueue.empty()){ 
    seed = seedQueue.front();
    seedQueue.pop_front();

    fileNum = seedQueue.front();
    seedQueue.pop_front();
    
    param = paramQueue.front();
    paramQueue.pop_front();
   
    MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, NEEDWORK,MPI_COMM_WORLD,&status);

    MPI_Send(&seed, 1, MPI_INT, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD); 
    MPI_Send(&fileNum, 1, MPI_INT,rank, 1, MPI_COMM_WORLD);
    MPI_Send(&param, 1, MPI_DOUBLE, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD); 
  }

  // Kill all other processors once all simulations are over
  for(rank = 1; rank < ntasks; rank++ ){
    MPI_Send(0,0,MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
    MPI_Send(0,0,MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
    MPI_Send(0,0,MPI_DOUBLE, rank, DIETAG, MPI_COMM_WORLD);
  } 


  
    
   return;
  } //end of master function


void slave(int rank){

   int seed, fileNum;
   double param;
   MPI_Status status;

   for(;;){

     // Receive RandSeed for simulation
     MPI_Recv(&seed,1,MPI_INT,MASTER,MPI_ANY_TAG,MPI_COMM_WORLD,&status); 
     MPI_Recv(&fileNum,1,MPI_INT,MASTER,MPI_ANY_TAG,MPI_COMM_WORLD,&status); 
     MPI_Recv(&param,1,MPI_DOUBLE,MASTER,MPI_ANY_TAG,MPI_COMM_WORLD,&status);

     if(status.MPI_TAG == DIETAG){
       return;
     }
  
     MPI_File file = files[fileNum]; 
     
     simulation(seed, rank, file, param);

     // Send back to Master that simulation is over
     MPI_Send(0,0,MPI_INT,MASTER,NEEDWORK,MPI_COMM_WORLD);
   }
}


void getSeeds(int& start, int& end, string& paramName, double& paramStart, double& paramStop, int& paramNum){
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

  // get randSeedStart
  while (strcmp(variable, "randSeedStart")) {
    input.getline(line, 256);
    variable = strtok(line, " ");
  }
  start = atoi(strtok(NULL, " "));

  // get randSeedEnd
  input.getline(line, 256);
  strtok(line, " ");
  end = atoi(strtok(NULL, " ")); 

  // get paramName
  while (strcmp(variable, "paramName")) {
    input.getline(line, 256);
    variable = strtok(line, " ");
  }
  paramName = strtok(NULL, " ");

  // get paramStart
  input.getline(line, 256);
  strtok(line, " ");
  paramStart = atof(strtok(NULL, " ")); 

  // get paramStop
  input.getline(line, 256);
  strtok(line, " ");
  paramStop = atof(strtok(NULL, " ")); 

  // get paramNum
  input.getline(line, 256);
  strtok(line, " ");
  paramNum = atoi(strtok(NULL, " ")); 
}

void openFiles(double paramStart, double paramStop, int paramNum) {

  double step = (paramStop - paramStart)/(paramNum - 1);
  double param = paramStart;
  for (int i = 0; i < paramNum; i++, param += step){
    
    MPI_File file;
    string file_name = "data1-" + to_string(param) + ".txt";
  
    MPI_File_open(MPI_COMM_WORLD, file_name.c_str(), MPI_MODE_WRONLY|MPI_MODE_CREATE, 
		MPI_INFO_NULL, &file);

    files[i] = file;

  } 
}

void closeFiles(int paramNum){
  
  for (int i = 0; i < paramNum; i++)
    MPI_File_close(&(files[i]));

}
