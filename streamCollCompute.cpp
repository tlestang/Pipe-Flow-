#include <iostream>
#ifndef __global__
#define __global__
#include "global.h"
#endif

void streamingAndCollisionComputeMacroBodyForce(double *fin, double *fout, double *rho, double *ux, double *uy, double beta, double tau)
  
{
  double rhoDum, ftemp;
  double uxDum, uyDum; double omega = 1./tau; double eu, eueu, u2, feq;
  double omega1 = 1.0-omega; double coeff_forcing = 1.0- 0.5*omega; double force_driving;
  int nx, ny;

    for(int x=0;x<Dx;x++)
    {
      for(int y=0;y<Dy;y++)
	{
	  rhoDum = 0.0;
	  uxDum = uyDum = 0.0;
	  for (int k=0;k<9;k++)
	    {
	      ftemp = fin[IDX(x,y,k)];
	      rhoDum += ftemp;
	      uxDum += ftemp*c[k][0];
	      uyDum += ftemp*c[k][1];
	    }
	  uxDum = uxDum/rhoDum + /*Influence of the force*/0.5*beta/rhoDum; uyDum /= rhoDum;
	  u2 = -1.5*(uxDum*uxDum + uyDum*uyDum);
	  for (int k=0;k<9;k++)
	    {
	      eu = c[k][0]*uxDum + c[k][1]*uyDum;
	      eueu = 4.5*eu*eu;
	      feq = w[k]*rhoDum*(1.0+3.0*eu+eueu+u2);

	      force_driving = w[k]*coeff_forcing*beta*(3.*(c[k][0]-uxDum)+9.*c[k][0]*eu);
	      fin[IDX(x,y,k)]= fin[IDX(x,y,k)]*omega1+feq*omega+force_driving;	  
	      /*Streaming*/
	      nx = (x + c[k][0]+Dx)%Dx; ny = (y + c[k][1]+Dy)%Dy;
	      fout[IDX(nx,ny,k)] = fin[IDX(x,y,k)];
	    }

	  /*Assign fields values to lattice arrays (heap)*/
	  ux[idx(x,y)] = uxDum; uy[idx(x,y)] = uyDum;
	  rho[idx(x,y)] = rhoDum;
	}
    }
}

void streamingAndCollisionComputeMacroBodyForceSpatial(double *fin, double *fout, double *rho, double *ux, double *uy, double beta0, double *map, double tau)
  
{
  double rhoDum, ftemp;
  double uxDum, uyDum; double omega = 1./tau; double eu, eueu, u2, feq;
  double omega1 = 1.0-omega; double coeff_forcing = 1.0- 0.5*omega; double force_driving;
  double beta;
  int nx, ny;

  double eps = 0.1;
    for(int x=0;x<Dx;x++)
    {
      for(int y=0;y<Dy;y++)
	{
	  beta = (1.0+eps*map[idx(x,y)])*beta0;
	  rhoDum = 0.0;
	  uxDum = uyDum = 0.0;
	  for (int k=0;k<9;k++)
	    {
	      ftemp = fin[IDX(x,y,k)];
	      rhoDum += ftemp;
	      uxDum += ftemp*c[k][0];
	      uyDum += ftemp*c[k][1];
	    }
	  uxDum = uxDum/rhoDum + /*Influence of the force*/0.5*beta/rhoDum; uyDum /= rhoDum;
	  u2 = -1.5*(uxDum*uxDum + uyDum*uyDum);
	  for (int k=0;k<9;k++)
	    {
	      eu = c[k][0]*uxDum + c[k][1]*uyDum;
	      eueu = 4.5*eu*eu;
	      feq = w[k]*rhoDum*(1.0+3.0*eu+eueu+u2);

	      force_driving = w[k]*coeff_forcing*beta*(3.*(c[k][0]-uxDum)+9.*c[k][0]*eu);
	      fin[IDX(x,y,k)]= fin[IDX(x,y,k)]*omega1+feq*omega+force_driving;	  
	      /*Streaming*/
	      nx = (x + c[k][0]+Dx)%Dx; ny = (y + c[k][1]+Dy)%Dy;
	      fout[IDX(nx,ny,k)] = fin[IDX(x,y,k)];
	    }

	  /*Assign fields values to lattice arrays (heap)*/
	  ux[idx(x,y)] = uxDum; uy[idx(x,y)] = uyDum;
	  rho[idx(x,y)] = rhoDum;
	}
    }
}



