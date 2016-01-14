#ifndef __global__
#define __global__
#include "global.h"
#endif
#include <stdlib.h>

void initializePopulations(double *fin, int Dx, int Dy)
{

  for(int x=0;x<Dx;x++)
    {
      for(int y=0;y<Dy;y++)
	{
	  for(int k=0;k<9;k++)
	    {
	      fin[IDX(x,y,k)] = w[k];
	    }
	}
    }
}

void initializeFields(double *fin, double *rho, double *ux, double *uy, int Dx, int Dy)
{
  double rho_, u, v, f;
  for(int x=0;x<Dx;x++)
    {
      for(int y=0;y<Dy;y++)
	{
	  u = v = rho_ = 0.0;
	  for(int k=0;k<9;k++)
	    {
	      f = fin[IDX(x,y,k)];
	      rho_ += f;
	      u += f*c[k][0];
	      v += f*c[k][1];
	    }
	  rho[idx(x,y)] = rho_;
	  ux[idx(x,y)] = u/rho_;
	  uy[idx(x,y)] = v/rho_;
      	}
    }
}



/*void initializePopulations(FILE* ifile, double ***fin, int Dx, int Dy)
{
  for(int x=0;x<Dx;x++)
    {
      for (int y=0;y<Dy;y++)
	{
	  for (int k=0;k<9;k++)
	    {
	      ifile >> fin[x][y][k];
	    }
	}
    }
    }*/

void initializeFields(double **rho, double ***u, int Dx, int Dy)
{
  double rho0 = 1.0; double zeroVelocity[2] = {0.0, 0.0};
  for(int x=0;x<Dx;x++)
    {
      for(int y=0;y<Dy;y++)
	{
	  rho[x][y] = rho0;
	  u[x][y][0] = zeroVelocity[0];
	  u[x][y][1] = zeroVelocity[1];
      	}
    }
}
