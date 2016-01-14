#ifndef __global__
#define __global__
#include "global.h"
#endif

void computeDomainNoSlipWalls_BB(double *fout, double *fin)
{
  /*Apply halfway BB to north and south wall, including corners*/
  int y0;
  for(int x=0;x<Dx;x++)
    {
      /*North boundary*/
      fout[IDX(x,Dy-1,4)] = fin[IDX(x,Dy-1,2)];
      fout[IDX(x,Dy-1,8)] = fin[IDX(x,Dy-1,6)];
      fout[IDX(x,Dy-1,7)] = fin[IDX(x,Dy-1,5)];

      /*South boundary*/
      fout[IDX(x,0,2)] = fin[IDX(x,0,4)];
      fout[IDX(x,0,5)] = fin[IDX(x,0,7)];
      fout[IDX(x,0,6)] = fin[IDX(x,0,8)];
    }
}


// void computeDomainNoSlipWalls_VERSA(double ***f, int Dx, int Dy)
// {
// /*This function apply Joris Verschaeve's no slip BC to the north and south wall*/
// /*NON OPTIMIZED VERSION*/
// /*UNUSABLE AT THE MOMENT BECAUSE IT USES OBSOLETE FUNCTION FEQ()*/
//   double fneq[9];
//   double uZero = 0.0;
//   int x0, y0;
//   double m_ps, rho;

//   /*North wall*/
//   y0 = Dy-1;
//   for(int x=1;x<Dx-1;x++)
//     {
//       m_ps = f[x][y0][6]+f[x][y0][2]+f[x][y0][5];
//       rho = 6.0*m_ps;
//       fneq[6] = f[x][y0][6] - feq(6,rho,uZero, uZero);
//       fneq[2] = f[x][y0][2] - feq(2,rho,uZero, uZero);
//       fneq[5] = f[x][y0][5] - feq(5,rho,uZero, uZero);
//       f[x][y0][0] = feq(0,rho,uZero, uZero); fneq[0] = 0.0;
//       f[x][y0][3] = feq(3,rho,uZero, uZero); fneq[3] = 0.0;
//       f[x][y0][1] = feq(1,rho,uZero, uZero); fneq[1] = 0.0;
//       fneq[7] = - fneq[6]; fneq[4] = - fneq[2]; fneq[8] = - fneq[5];
//       f[x][y0][7] = feq(7,rho,uZero, uZero) + fneq[7]; 
//       f[x][y0][4] = feq(4,rho,uZero, uZero) + fneq[4]; 
//       f[x][y0][8] = feq(8,rho,uZero, uZero) + fneq[8];
//     }

//   /*South wall*/
//   y0 = 0;
//   for(int x=1;x<Dx-1;x++)
//     {			   
//       m_ps = f[x][y0][7]+f[x][y0][8]+f[x][y0][4];
//       rho = 6.0*m_ps;
//       fneq[7] = f[x][y0][7] - feq(7,rho,uZero, uZero);
//       fneq[4] = f[x][y0][4] - feq(4,rho,uZero, uZero);
//       fneq[8] = f[x][y0][8] - feq(8,rho,uZero, uZero);
//       f[x][y0][0] = feq(0,rho,uZero, uZero); fneq[0] = 0.0;
//       f[x][y0][3] = feq(3,rho,uZero, uZero); fneq[3] = 0.0;
//       f[x][y0][1] = feq(1,rho,uZero, uZero); fneq[1] = 0.0;
//       fneq[6] = - fneq[7]; fneq[2] = - fneq[4]; fneq[5] = - fneq[8];
//       f[x][y0][6] = feq(6,rho,uZero, uZero) + fneq[6]; 
//       f[x][y0][2] = feq(2,rho,uZero, uZero) + fneq[2]; 
//       f[x][y0][5] = feq(5,rho,uZero, uZero) + fneq[5];
//     }
  
//   /*Top left corner*/
//   x0=0; y0=Dy-1;
//   m_ps=f[x0][y0][6];
//   rho = 36.0*m_ps;
//   for(int k=0;k<9;k++){f[x0][y0][k]=feq(k,rho,uZero, uZero);}

//   /*Top right corner*/
//   x0=Dx-1; y0=Dy-1;
//   m_ps=f[x0][y0][5];
//   rho = 36.0*m_ps;
//   for(int k=0;k<9;k++){f[x0][y0][k]=feq(k,rho,uZero, uZero);}

//   /*Bottom right corner*/
//   x0=Dx-1; y0=0;
//   m_ps=f[x0][y0][8];
//   rho = 36.0*m_ps;
//   for(int k=0;k<9;k++){f[x0][y0][k]=feq(k,rho,uZero, uZero);}

//   /*Bottom left corner*/
//   x0=0; y0=0;
//   m_ps=f[x0][y0][7];
//   rho = 36.0*m_ps;
//   for(int k=0;k<9;k++){f[x0][y0][k]=feq(k,rho,uZero, uZero);}
// }
