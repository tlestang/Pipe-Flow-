#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <malloc.h>
#include <unistd.h>
#include "stdio.h"
#include <sys/time.h>

#ifndef __global__
#define __global__
#include "global.h"
#endif

#include "initialize_lattice_arrays.h"
#include "streamCollCompute.h"
#include "boundaryConditions.h"
#include "force.h"
#include "write_vtk.h"

using namespace std;

int c[9][2] = {{0,0}, {1,0}, {0,1}, {-1,0}, {0,-1}, {1,1}, {-1,1}, {-1,-1}, {1,-1}};
double w[9]={4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
int Dx, Dy, xmin, xmax, ymin, ymax;

int main()
{
  double ot = 1./3;
  double entry;
  /*Parameters for LB simulation*/
  int nbOfChunks, nbOfTimeSteps, numberOfTransientSteps, Lx, Ly;
  int facquVtk, facquRe, facquForce, facquU;
  double tau, beta;
  double Ma;   //Mach number
  string folderName, inputPopsFileName;
  /*Reads input file*/
  ifstream input_file("input.datin");
  input_file >> nbOfChunks;
  input_file >> nbOfTimeSteps;
  input_file >> Lx; Ly = Lx;
  input_file >> tau;
  input_file >> Ma;
  input_file >> folderName;
  input_file >> inputPopsFileName;
  input_file >> facquVtk;
  input_file >> facquForce;
  input_file >> facquU;
  input_file.close();
  
  /*Compute or define other parameters*/
  Dy = 4*Ly + 1, Dx = Dy; // 2*(Dy-1) + 1;
  xmin = (Dx-1)/2; xmax = xmin + Lx;
  ymin = (Dy-1)/2 - Ly/2; ymax = ymin + Ly;
  double cs = 1./sqrt(3); double rho0 = 1.0;
  double u0 = cs*cs*Ma;
  double nu = ot*(tau-0.5);
  double omega = 1.0/tau;
  beta = 8*nu*u0/((Dy-1)/2)/((Dy-1)/2);
  

  //----------- Misc ----------
  double uxSum = 0.0, uxMean;
  double F; int tt=0;
  int dummy, dummy2;
  /*Populations and macroscopic fields*/
  double *fin, *fout, *rho, *ux, *uy, *temp;
  /* --- | Create folder for storing data | ---  */
  string instru = "mkdir " + folderName;
  system(instru.c_str());
  instru = "mkdir " + folderName + "/vtk_fluid/";
  system(instru.c_str());


  /* --- | Create parameters file | --- */
  string openParamFile = folderName + "/parameters.datout";
  ofstream param;
  param.open(openParamFile.c_str());
  param << "Number of timesteps : " << nbOfChunks << "X" << nbOfTimeSteps << endl;
  param << "L : "  << Lx << endl;
  param << "Dx : " << Dx << endl;
  param << "Dy : " << Dy << endl;
  param << "tau : " << tau << endl;
  param << "beta : " << beta << endl;
  param.close();

  /* ---- | Allocate populations and fields | --- */

  fin = (double *) memalign(getpagesize(), Dx*Dy*9*sizeof(double));
  fout = (double *) memalign(getpagesize(), Dx*Dy*9*sizeof(double));
  rho = (double *) memalign(getpagesize(), Dx*Dy*sizeof(double));
  ux = (double *) memalign(getpagesize(), Dx*Dy*sizeof(double));
  uy = (double *) memalign(getpagesize(), Dx*Dy*sizeof(double));
  
  if(inputPopsFileName != "0")
    {
      ifstream popFile(inputPopsFileName.c_str(), ios::binary);
      cout << "Initialized populations taken from " << inputPopsFileName << endl;
      popFile.read((char*)&fin[0], Dx*Dy*9*sizeof(double));
      popFile.close();
    }
  else
    {
  /*Initialization of population to equilibrium value*/
      cout << "Initializing pops to equilibrium value" << endl;
      initializePopulations(fin, Dx, Dy);
      initializeFields(fin, rho, ux, uy, Dx, Dy);
    }
  
  /*Initialize counters*/
  dummy = 0; dummy2 = 0;

  /*Open output files for Reynolds, Mach and force*/

  string openForceFile = folderName + "/data_force.datout";
  string openuxFile = folderName + "/ux_t.datout";
  ofstream ReFile, MaFile, forceFile, uxFile;

  uxFile.open(openuxFile.c_str(), ios::binary);
  forceFile.open(openForceFile.c_str(), ios::binary);

  /*Start LBM*/

  // for(int chunkID=0;chunkID<nbOfChunks;chunkID++)
  // {
  //        if(chunkID%(nbOfChunks/100)==0){dummy2++; cout<<"Running : " << dummy2<<"%"<<endl;/*\r"; fflush(stdout);*/}

  for (int lbTimeStepCount=0; lbTimeStepCount<nbOfTimeSteps;lbTimeStepCount++)
    {
      if(lbTimeStepCount%(nbOfTimeSteps/100)==0)
       	dummy2++; cout<<dummy2<<"%\r"; fflush(stdout);
      if(lbTimeStepCount%facquVtk==0)
      	{
      	  write_fluid_vtk(tt, Dx, Dy, rho, ux, uy, folderName.c_str());
      	  tt++;
      	}
      /*Collision and streaming - Macroscopic fields*/
      //streamingAndCollisionComputeMacroBodyForce(popHeapIn, popHeapOut, rhoHeap, uFieldHeap, Dx, Dy, tau, beta);
      streamingAndCollisionComputeMacroBodyForce(fin, fout, rho, ux, uy, beta, tau);
      computeDomainNoSlipWalls_BB(fout, fin);
      computeSquareBounceBack_TEST(fout, fin);
      /*Reset square nodes to equilibrium*/
      for(int x=xmin+1;x<xmax;x++)
	{
	  for(int y=ymin+1;y<ymax;y++)
	    {
	      for(int k=0;k<9;k++)
		{
		  fout[IDX(x,y,k)] = w[k];
		}
	    }
	}
      /*Swap populations*/
      temp = fin;
      fin = fout;
      fout = temp;

      /*Compute and Write force on disk*/
      if(lbTimeStepCount%facquForce==0)
	{
	  F = computeForceOnSquare(fin, omega);
	  forceFile.write((char*)&F, sizeof(double));
	}


      if(lbTimeStepCount%facquU==0)
      	{      
      	  uxFile.write((char*)&ux[idx(Dx/4,Dy/4)], sizeof(double));
      	}
      
    }
  //}
 
  forceFile.close();
  uxFile.close();
  /*End of run - Save populations on disk*/
  /*and complete parameters file*/
  string popsFileName = folderName + "/pops.datout";
  ofstream pops_output_file(popsFileName.c_str(), ios::binary);
  pops_output_file.write((char*)&fin[0], Dx*Dy*9*sizeof(double));
  pops_output_file.close();
}
