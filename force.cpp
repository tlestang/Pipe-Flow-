#include <iostream>

#ifndef __global__
#define __global__
#include "global.h"
#endif

double computeForceOnSquare(double *f, double omega)
{
  double ot = 1./3.;
  int x0, y0;
  double eu, eueu, u2;
  double fEast, fWest, fNorth, fSouth;
  double rho_, ux, uy, Pi_xx, Pi_xy;
  double fneq, ftemp, feq;
  double totalForce;
  double coeff_force = 1.0-.5*omega;
  
  /*West side*/
  /*Compute Pi_xx*/
  x0 = xmin;
  fWest = 0.0;
  for(int y=ymin;y<ymax+1;y++)
    {
      /*Compute local macro. fields*/
      rho_ = 0.0; ux = 0.0; uy = 0.0;
      for(int k=0;k<9;k++)
	{
	  ftemp = f[IDX(x0,y,k)];
	  rho_ += ftemp;
	  ux += ftemp*c[k][0];
	  uy += ftemp*c[k][1];
	}
      uy /= rho_;
      ux /= rho_;
      /*Compute tensor Pi1*/
      Pi_xx = 0.0;
      u2 = -1.5*(ux*ux + uy*uy);
      for(int k=0;k<9;k++)
	{
	  eu = c[k][0]*ux + c[k][1]*uy;
	  eueu = 4.5*eu*eu;
	  feq = w[k]*rho_*(1.0+3.0*eu+eueu+u2);
	  fneq = f[IDX(x0,y,k)] - feq;
	  Pi_xx += fneq*c[k][0]*c[k][0];
	}
      /*Compute force*/
      fWest += rho_*ot + coeff_force*Pi_xx;
    }

  /*East side*/
  x0 = xmax;
  fEast = 0.0;
  for(int y=ymin;y<ymax+1;y++)
    {
      /*Compute local macro. fields*/
      rho_ = 0.0; ux = 0.0; uy = 0.0;
      for(int k=0;k<9;k++)
	{
	  ftemp = f[IDX(x0,y,k)];
	  rho_ += ftemp;
	  ux += ftemp*c[k][0];
	  uy += ftemp*c[k][1];
	}
      uy /= rho_;
      ux /= rho_;
      /*Compute tensor Pi1*/
      Pi_xx = 0.0;
      u2 = -1.5*(ux*ux + uy*uy);
      for(int k=0;k<9;k++)
	{
	  eu = c[k][0]*ux + c[k][1]*uy;
	  eueu = 4.5*eu*eu;
	  feq = w[k]*rho_*(1.0+3.0*eu+eueu+u2);
	  fneq = f[IDX(x0,y,k)] - feq;
	  Pi_xx += fneq*c[k][0]*c[k][0];
	}
      fEast += - rho_*ot - coeff_force*Pi_xx;
    }

  /*North side*/
  y0 = ymax;
  fNorth = 0.0;
  for(int x=xmin;x<xmax+1;x++)
    {
      /*Compute local macro. fields*/
      rho_ = 0.0; ux = 0.0; uy = 0.0;
      for(int k=0;k<9;k++)
	{
	  ftemp = f[IDX(x,y0,k)];
	  rho_ += ftemp;
	  ux += ftemp*c[k][0];
	  uy += ftemp*c[k][1];
	}
      uy /= rho_;
      ux /= rho_;
      /*Compute tensor Pi1*/
      Pi_xy = 0.0;
      u2 = -1.5*(ux*ux + uy*uy);
      for(int k=0;k<9;k++)
	{
	  eu = c[k][0]*ux + c[k][1]*uy;
	  eueu = 4.5*eu*eu;
	  feq = w[k]*rho_*(1.0+3.0*eu+eueu+u2);
	  fneq = f[IDX(x,y0,k)] - feq;
	  Pi_xy += fneq*c[k][0]*c[k][1];
	}
      fNorth += - coeff_force*Pi_xy;
    }

  /*South side*/
  y0 = ymin;
  fSouth = 0.0;
  for(int x=xmin;x<xmax+1;x++)
    {
      rho_ = 0.0; ux = 0.0; uy = 0.0;
      /*Compute local macro. fields*/
      for(int k=0;k<9;k++)
	{
	  ftemp = f[IDX(x,y0,k)];
	  rho_ += ftemp;
	  ux += ftemp*c[k][0];
	  uy += ftemp*c[k][1];
	}
      uy /= rho_;
      ux /= rho_;
      /*Compute tensor Pi1*/
      Pi_xy = 0.0;
      u2 = -1.5*(ux*ux + uy*uy);
      for(int k=0;k<9;k++)
	{
	  eu = c[k][0]*ux + c[k][1]*uy;
	  eueu = 4.5*eu*eu;
	  feq = w[k]*rho_*(1.0+3.0*eu+eueu+u2);
	  fneq = f[IDX(x,y0,k)] - feq;
	  Pi_xy += fneq*c[k][0]*c[k][1];
	}
      fSouth += + coeff_force*Pi_xy;
    }

  totalForce = fWest + fEast + fNorth + fSouth;
  return totalForce;

}


