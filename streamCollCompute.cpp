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
	      uxDum += ftemp*e[k][0];
	      uyDum += ftemp*e[k][1];
	    }
	  uxDum = uxDum/rhoDum + /*Influence of the force*/0.5*beta/rhoDum; uyDum /= rhoDum;
	  u2 = -1.5*(uxDum*uxDum + uyDum*uyDum);
	  for (int k=0;k<9;k++)
	    {
	      eu = e[k][0]*uxDum + e[k][1]*uyDum;
	      eueu = 4.5*eu*eu;
	      feq = w[k]*rhoDum*(1.0+3.0*eu+eueu+u2);
	      force_driving = w[k]*coeff_forcing*beta*(3.*(e[k][0]-uxDum)+9.*e[k][0]*eu);
	      fin[IDX(x,y,k)]= fin[IDX(x,y,k)]*omega1+feq*omega+force_driving;	  
	      /*Streaming*/
	      nx = (x + e[k][0]+Dx)%Dx; ny = (y + e[k][1]+Dy)%Dy;
	      fout[IDX(nx,ny,k)] = fin[IDX(x,y,k)];
	    }

	  /*Assign fields values to lattice arrays (heap)*/
	  ux[idx(x,y)] = uxDum; uy[idx(x,y)] = uyDum;
	  rho[idx(x,y)] = rhoDum;
	}
    }
}

void streamingAndCollisionComputeMacro
(double ***fin, double ***fout, double **rho, double ***u, int Dx, int Dy, double tau)
  
{
  double w[9]={4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
  double rhoDum, ftemp;
  double ux, uy; double omega = 1./tau; double eu, eueu, u2, feq;
  double omega1 = 1.0-omega; 
  int e[9][2] = {{0,0}, {1,0}, {0,1}, {-1,0}, {0,-1}, {1,1}, {-1,1}, {-1,-1}, {1,-1}};
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
	      ux += ftemp*e[k][0];
	      uy += ftemp*e[k][1];
	    }
	  ux = ux/rhoDum; uy /= rhoDum;
	  u2 = -1.5*(ux*ux + uy*uy);
	  for (int k=0;k<9;k++)
	    {
	      eu = e[k][0]*ux + e[k][1]*uy;
	      eueu = 4.5*eu*eu;
	      feq = w[k]*rhoDum*(1.0+3.0*eu+eueu+u2);
	      fin[x][y][k] = fin[x][y][k]*omega1+feq*omega;	  
	      /*Streaming*/
	      nx = (x + e[k][0]+Dx)%Dx; ny = (y + e[k][1]+Dy)%Dy;
	      fout[nx][ny][k] = fin[x][y][k];
	    }

	  /*Assign fields values to lattice arrays (heap)*/
	  u[x][y][0] = ux; u[x][y][1] = uy;
	  rho[x][y] = rhoDum;
	}
    }
}

/*void streamingAndCollisionComputeSquare
(double ***fin, double ***fout, double **rho, double ***u, int Dx, int Dy, double tau, int Lx, int Ly)
  
{
  double rhoDum;
  double ux, uy; double omega = 1./tau; 
  int xmin = (Dx-1)/3; int xmax = xmin+Lx-1;
  int ymin = Dy/2 - Ly/2; int ymax = ymin + Ly -1;
  int e[9][2] = {{0,0}, {1,0}, {0,1}, {-1,0}, {0,-1}, {1,1}, {-1,1}, {-1,-1}, {1,-1}};
  int nx, ny;
  
  for(int x=0;x<Dx;x++)
    {
      for(int y=0;y<Dy;y++)
	{
	  if(x<=xmin || x>=xmax || y<=ymin || y>= ymax)
	      {
	  rhoDum = 0.0;
	  ux = uy = 0.0;
	  for (int k=0;k<9;k++)
	    {
	      rhoDum += fin[x][y][k];
	      ux += fin[x][y][k]*e[k][0];
	      uy += fin[x][y][k]*e[k][1];
	    }
	  ux /= rhoDum; uy /= rhoDum;

	  for (int k=0;k<9;k++)
	    {
	      // Collision
	  fin[x][y][k] *= (1.0-omega);
	  fin[x][y][k] += omega*feq(k, rhoDum, ux, uy);
	  //Streaming
	  nx = (x + e[k][0]+Dx)%Dx; ny = (y + e[k][1]+Dy)%Dy;
	  fout[nx][ny][k] = fin[x][y][k];
	    }
	  //Assign fields values to lattice arrays (heap)
	  u[x][y][0] = ux; u[x][y][1] = uy;
	  rho[x][y] = rhoDum;
	      }
	}
    }
    }*/

//void streamingAndCollisionComputeSquare
//(double ***fin, double ***fout, double **rho, double ***u, int Dx, int Dy, double tau)
/*This function applies the streaming and collision operator to the fluid AROUND the square. It considers the inside the square to be NONFLUID*/
/*This function is UNUSABLE at the moment because it uses the OBSOLETE function feq()*/
/*This function is NON OPTIMIZED*/
//{
  // double rhoDum;
  // double ux, uy; double omega = 1./tau;
  // int e[9][2] = {{0,0}, {1,0}, {0,1}, {-1,0}, {0,-1}, {1,1}, {-1,1}, {-1,-1}, {1,-1}};
  // int nx, ny;

//    for(int x=0;x<Dx;x++)
//     {
//       for(int y=0;y<Dy;y++)
// 	{
// 	  rhoDum = 0.0;
// 	  ux = uy = 0.0;
// 	  for (int k=0;k<9;k++)
// 	    {
// 	      rhoDum += fin[x][y][k];
// 	      ux += fin[x][y][k]*e[k][0];
// 	      uy += fin[x][y][k]*e[k][1];
// 	    }
// 	  ux /= rhoDum; uy /= rhoDum;

// 	  for (int k=0;k<9;k++)
// 	    {
// 	  /* Collision*/
// 	  fin[x][y][k] *= (1.0-omega);
// 	  fin[x][y][k] += omega*feq(k, rhoDum, ux, uy);
// 	  /*Streaming*/
// 	  nx = (x + e[k][0]+Dx)%Dx; ny = (y + e[k][1]+Dy)%Dy;
// 	  fout[nx][ny][k] = fin[x][y][k];
// 	    }
// 	  /*Assign fields values to lattice arrays (heap)*/
// 	  u[x][y][0] = ux; u[x][y][1] = uy;
// 	  rho[x][y] = rhoDum;
// 	}
//     }
// }*/



