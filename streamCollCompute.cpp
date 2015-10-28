#include "equilibrium.h"
#include <iostream>

void streamingAndCollisionComputeSquareBodyForce
(double ***fin, double ***fout, double **rho, double ***u, int Dx, int Dy, double tau, int xmin, int xmax, int ymin, int ymax, double beta)
  
{
  double w[9]={4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
  double rhoDum, ftemp;
  double ux, uy; double omega = 1./tau; double eu, eueu, u2, feq;
  double omega1 = 1.0-omega; double coeff_forcing = 1.0- 0.5*omega; double force_driving;
  int c[9][2] = {{0,0}, {1,0}, {0,1}, {-1,0}, {0,-1}, {1,1}, {-1,1}, {-1,-1}, {1,-1}};
  int nx, ny;
  
  for(int x=0;x<Dx;x++)
    {
      for(int y=0;y<Dy;y++)
	{
	  rhoDum = 0.0;
	  ux = uy = 0.0;
	  for (int k=0;k<9;k++)
	    {
	      ftemp = fin[x][y][k];
	      rhoDum += ftemp;
	      ux += ftemp*c[k][0];
	      uy += ftemp*c[k][1];
	    }
	  ux = ux/rhoDum + /*Influence of the force*/0.5*beta; uy /= rhoDum;
	  u2 = -1.5*(ux*ux + uy*uy);
	  for (int k=0;k<9;k++)
	    {
	      eu = c[k][0]*ux + c[k][1]*uy;
	      eueu = 4.5*eu*eu;
	      feq = w[k]*rhoDum*(1.0+3.0*eu+eueu+u2);
	      force_driving = w[k]*coeff_forcing*beta*(3.*(c[k][0]-ux)+9.*c[k][0]*eu);
	      fin[x][y][k] = fin[x][y][k]*omega1+feq*omega+force_driving;	  
	      /*Streaming*/
	      nx = (x + c[k][0]+Dx)%Dx; ny = (y + c[k][1]+Dy)%Dy;
	      fout[nx][ny][k] = fin[x][y][k];
	    }

	  /*Assign fields values to lattice arrays (heap)*/
	  u[x][y][0] = ux; u[x][y][1] = uy;
	  rho[x][y] = rhoDum;
	}
    }
}

