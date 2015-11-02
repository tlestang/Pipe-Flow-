#include <iostream>
using namespace std;

double computeForceOnSquare(double ***f, int xmax, int xmin, int ymax, int ymin, double omega)
{

  double ot = 1./3.;
  double w[9]={4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
  int e[9][2] = {{0,0}, {1,0}, {0,1}, {-1,0}, {0,-1}, {1,1}, {-1,1}, {-1,-1}, {1,-1}};
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
	  ftemp = f[x0][y][k];
	  rho_ += ftemp;
	  ux += ftemp*e[k][0];
	  uy += ftemp*e[k][1];
	}
      uy /= rho_;
      ux /= rho_;
      /*Compute tensor Pi1*/
      Pi_xx = 0.0;
      u2 = -1.5*(ux*ux + uy*uy);
      for(int k=0;k<9;k++)
	{
	  eu = e[k][0]*ux + e[k][1]*uy;
	  eueu = 4.5*eu*eu;
	  feq = w[k]*rho_*(1.0+3.0*eu+eueu+u2);
	  fneq = f[x0][y][k] - feq;
	  Pi_xx += fneq*e[k][0]*e[k][0];
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
	  ftemp = f[x0][y][k];
	  rho_ += ftemp;
	  ux += ftemp*e[k][0];
	  uy += ftemp*e[k][1];
	}
      uy /= rho_;
      ux /= rho_;
      /*Compute tensor Pi1*/
      Pi_xx = 0.0;
      u2 = -1.5*(ux*ux + uy*uy);
      for(int k=0;k<9;k++)
	{
	  eu = e[k][0]*ux + e[k][1]*uy;
	  eueu = 4.5*eu*eu;
	  feq = w[k]*rho_*(1.0+3.0*eu+eueu+u2);
	  fneq = f[x0][y][k] - feq;
	  Pi_xx += fneq*e[k][0]*e[k][0];
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
	  ftemp = f[x][y0][k];
	  rho_ += ftemp;
	  ux += ftemp*e[k][0];
	  uy += ftemp*e[k][1];
	}
      uy /= rho_;
      ux /= rho_;
      /*Compute tensor Pi1*/
      Pi_xy = 0.0;
      u2 = -1.5*(ux*ux + uy*uy);
      for(int k=0;k<9;k++)
	{
	  eu = e[k][0]*ux + e[k][1]*uy;
	  eueu = 4.5*eu*eu;
	  feq = w[k]*rho_*(1.0+3.0*eu+eueu+u2);
	  fneq = f[x][y0][k] - feq;
	  Pi_xy += fneq*e[k][0]*e[k][1];
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
	  ftemp = f[x][y0][k];
	  rho_ += ftemp;
	  ux += ftemp*e[k][0];
	  uy += ftemp*e[k][1];
	}
      uy /= rho_;
      ux /= rho_;
      /*Compute tensor Pi1*/
      Pi_xy = 0.0;
      u2 = -1.5*(ux*ux + uy*uy);
      for(int k=0;k<9;k++)
	{
	  eu = e[k][0]*ux + e[k][1]*uy;
	  eueu = 4.5*eu*eu;
	  feq = w[k]*rho_*(1.0+3.0*eu+eueu+u2);
	  fneq = f[x][y0][k] - feq;
	  Pi_xy += fneq*e[k][0]*e[k][1];
	}
      fSouth += + coeff_force*Pi_xy;
    }

  totalForce = fWest + fEast + fNorth + fSouth;

  return totalForce;

}


