#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <sys/time.h>

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
  double tau, beta;
  double Ma;   //Mach number
  /*Reads input file*/
  ifstream input_file("input.datin");
  input_file >> nbOfChunks;
  input_file >> nbOfTimeSteps;
  input_file >> Lx; Ly = Lx;
  input_file >> tau;
  input_file >> Ma;
  input_file.close();
  /*Compute or define other parameters*/
  int Dy = 4*Ly + 1, Dx = 2*(Dy-1) + 1;
  int xmin = (Dx-1)/2; int xmax = xmin + Lx;
  int ymin = (Dy-1)/2 - Ly/2; int ymax = ymin + Ly;
  double cs = 1./sqrt(3); double rho0 = 1.0;
  double u0 = cs*cs*Ma;
  double nu = ot*(tau-0.5);
  double omega = 1.0/tau;
  beta = 8*nu*u0/((Dy-1)/2)/((Dy-1)/2);
  

  //----------- Misc ----------
  double w[9]={4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
  double uxSum = 0.0, uxMean;
  double F; int tt=0;
  int dummy, dummy2;
  /*Populations and macroscopic fields*/
  double ***popHeapIn, ***popHeapOut, ***uFieldHeap, ***temp;
  double **rhoHeap;

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


  /*Initialization of population to equilibrium value*/
      cout << "Initializing pops to equilibrium value" << endl;
  initializePopulations(popHeapIn, Dx, Dy);
  initializeFields(rhoHeap, uFieldHeap, Dx, Dy);
  
  /*Initialize counters*/
  dummy = 0; dummy2 = 0;

//Variables for performance evaluation
struct timeval start, end;
                                  /* --- START LBM --- */

 //TICK TIMER
      gettimeofday(&start,NULL);
  for (int lbTimeStepCount=0; lbTimeStepCount<nbOfTimeSteps;lbTimeStepCount++)
    {
      
      
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
      
    }
  
  gettimeofday(&end,NULL);
  double t = (end.tv_sec - start.tv_sec)*1e6 + (end.tv_usec - start.tv_usec);
  cout << "Average time per ts in us : " << t/nbOfTimeSteps << endl;
  cout << "Total runtime in s : " << t/(1e6) << endl;
}
