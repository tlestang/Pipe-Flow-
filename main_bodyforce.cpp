#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "stdio.h"

#include "initialize_lattice_arrays.h"
#include "streamCollCompute.h"
#include "boundaryConditions.h"
#include "force.h"
#include "write_vtk.h"

using namespace std;

int main()
{
  double ot = 1./3;
  double entry;
  /*Parameters for LB simulation*/
  int nbOfChunks, nbOfTimeSteps, numberOfTransientSteps, Lx, Ly;
  int facquVtk, facquRe, facquForce;
  double tau, beta;
  string folderName, inputPopsFileName;
  /*Reads input file*/
  ifstream input_file("input.datin");
  input_file >> nbOfChunks;
  input_file >> nbOfTimeSteps;
  input_file >> Lx; Ly = Lx;
  input_file >> tau;
  input_file >> beta;
  input_file >> folderName;
  input_file >> inputPopsFileName;
  input_file >> facquVtk;
  input_file >> facquRe;
  input_file >> facquForce;
  input_file.close();
  /*Compute or define other parameters*/
  int Dy = 4*Ly, Dx = 2*Dy;
  int xmin = (Dx-1)/2; int xmax = xmin + Lx;
  int ymin = Dy/2 - Ly/2; int ymax = ymin + Ly;
  double cs = 1./sqrt(3); double rho0 = 1.0;

  double Ma;   //Mach number
  double nu = ot*(tau-0.5);
  double omega = 1.0/tau;
  //----------- Misc ----------
  double w[9]={4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
  double uxSum = 0.0, uxMean;
  double F; int tt=0;
  int dummy, dummy2;
  /*Populations and macroscopic fields*/
  double ***popHeapIn, ***popHeapOut, ***uFieldHeap, ***temp;
  double **rhoHeap;
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

  popHeapIn = new double**[Dx]; popHeapOut = new double**[Dx];
  for (int i=0;i<Dx;i++)
    {
      popHeapIn[i] = new double*[Dy];
      popHeapOut[i] = new double*[Dy];
      for(int j=0;j<Dy;j++)
	{
	  popHeapOut[i][j] = new double[9];
	  popHeapIn[i][j] = new double[9];
	}
    }
  rhoHeap = new double*[Dx]; uFieldHeap = new double**[Dx];
  for(int i=0;i<Dx;i++)
    {
      rhoHeap[i] = new double[Dy];
      uFieldHeap[i] = new double*[Dy];
      for (int j=0;j<Dy;j++)
	{
	  uFieldHeap[i][j] = new double[2];
	}
    }


  if(inputPopsFileName != "0")
    {
      ifstream popFile(inputPopsFileName.c_str());
      cout << "Initialized populations taken from " << inputPopsFileName << endl;
      for(int x=0;x<Dx;x++)
	{
	  for(int y=0;y<Dy;y++)
	    {
	      for(int k=0;k<9;k++)
		{
		  popFile >> popHeapIn[x][y][k];
		}
	    }
	}
      popFile.close();
    }
  else
    {
  /*Initialization of population to equilibrium value*/
      cout << "Initializing pops to equilibrium value" << endl;
  initializePopulations(popHeapIn, Dx, Dy);
  initializeFields(rhoHeap, uFieldHeap, Dx, Dy);
    }
  
  /*Initialize counters*/
  dummy = 0; dummy2 = 0;

  /*Open output files for Reynolds, Mach and force*/
  string openReFile = folderName + "/re_t.datout";
  string openMaFile = folderName + "/ma_t.datout";
  string openForceFile = folderName + "/data_force.datout";
  ofstream ReFile, MaFile, data_force;
  ReFile.open(openReFile.c_str());
  MaFile.open(openMaFile.c_str());
  data_force.open(openForceFile.c_str());

  /*Start LBM*/
for(int chunkID=0;chunkID<nbOfChunks;chunkID++)
{
  if(chunkID%(nbOfChunks/100)==0){dummy2++; cout<<"Running : " << dummy2<<"%"<<endl;/*\r"; fflush(stdout);*/}

    for (int lbTimeStepCount=0; lbTimeStepCount<nbOfTimeSteps;lbTimeStepCount++)
    {
      /*if(lbTimeStepCount%(nbOfTimeSteps/100)==0)
	dummy++; cout<<dummy<<"%\r"; fflush(stdout);*/
      if(lbTimeStepCount%facquVtk==0)
	{
	write_fluid_vtk(tt, Dx, Dy, rhoHeap, uFieldHeap, folderName.c_str());
	tt++;
	}
      /*Collision and streaming - Macroscopic fields*/
      streamingAndCollisionComputeMacroBodyForce(popHeapIn, popHeapOut, rhoHeap, uFieldHeap, Dx, Dy, tau, beta);
      computeDomainNoSlipWalls_BB(popHeapOut, popHeapIn, Dx, Dy);
      computeSquareBounceBack_TEST(popHeapOut, popHeapIn, xmin, xmax, ymin, ymax);
      /*Reset square nodes to equilibrium*/
      for(int x=xmin+1;x<xmax;x++)
	{
	  for(int y=ymin+1;y<ymax;y++)
	    {
	      for(int k=0;k<9;k++)
		{
		  popHeapOut[x][y][k] = w[k];
		}
	    }
	}
      /*Swap populations*/
      temp = popHeapIn;
      popHeapIn = popHeapOut;
      popHeapOut = temp;

      /*Compute and Write force on disk*/
      if(lbTimeStepCount%facquForce==0)
      {
	F = computeForceOnSquare(popHeapIn, xmax, xmin, ymax, ymin, omega);
	data_force << F << endl;
	}
      /*Compute Reynolds number*/
      if(lbTimeStepCount%facquRe==0)
	{      
	  for(int y=0;y<Dy;y++)
	    {
	      uxSum += uFieldHeap[Dx/4][y][0];
	    }
	  uxMean = uxSum/Dy;
	  ReFile << lbTimeStepCount + chunkID*nbOfTimeSteps << " " << (2.0*uxMean*(Ly-1))/nu << endl;
	  MaFile << lbTimeStepCount + chunkID*nbOfTimeSteps << " " << uxMean/cs << endl;
	  uxSum=0.0; 
	  }
    }
 }
 ReFile.close();
 MaFile.close();
 data_force.close();
 /*End of run - Save populations on disk*/
 /*and complete parameters file*/
 string popsFileName = folderName + "/pops.datout";
 ofstream pops_output_file(popsFileName.c_str());
  for(int x=0;x<Dx;x++)
   {
     for(int y=0;y<Dy;y++)
       {
	 for(int k=0;k<9;k++)
	   {
	     pops_output_file << popHeapIn[x][y][k] << endl;
	   }
       }
   }
   pops_output_file.close();
}
