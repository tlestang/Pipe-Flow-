#ifndef __global__
#define __global__
#include "global.h"
#endif

void computeSquareBounceBack_FULLSTREAM(double ***fout, double ***fin, int xmin, int xmax, int ymin, int ymax)
{
  int opp[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
  for(int k=0;k<9;k++)
    {
      for(int x=xmin;x<xmax+1;x++)
	{
	  fout[x][ymax][k] = fin[x][ymax][opp[k]];
	  fout[x][ymin][k] = fin[x][ymin][opp[k]];
	}
            for(int y=ymin;y<ymax+1;y++)
	{
	  fout[xmax][y][k] = fin[xmax][y][opp[k]];
	  fout[xmin][y][k] = fin[xmin][y][opp[k]];
	}
    }
}

void computeSquareBounceBack_TEST(double *fout, double *fin)
{
    for(int x=xmin+1;x<xmax;x++)
    {
      fout[IDX(x,ymin,7)] = fin[IDX(x,ymin,5)];
      fout[IDX(x,ymin,4)] = fin[IDX(x,ymin,2)];
      fout[IDX(x,ymin,8)] = fin[IDX(x,ymin,6)];

      fout[IDX(x,ymax,5)] = fin[IDX(x,ymax,7)];
      fout[IDX(x,ymax,2)] = fin[IDX(x,ymax,4)];
      fout[IDX(x,ymax,6)] = fin[IDX(x,ymax,8)];
    }
  for(int y=ymin+1;y<ymax;y++)
    {
      fout[IDX(xmin,y,6)] = fin[IDX(xmin,y,8)];
      fout[IDX(xmin,y,3)] = fin[IDX(xmin,y,1)];
      fout[IDX(xmin,y,7)] = fin[IDX(xmin,y,5)];

      fout[IDX(xmax,y,8)] = fin[IDX(xmax,y,6)];
      fout[IDX(xmax,y,1)] = fin[IDX(xmax,y,3)];
      fout[IDX(xmax,y,5)] = fin[IDX(xmax,y,7)];
    }
  fout[IDX(xmin,ymin,7)] = fin[IDX(xmin,ymin,5)];
  fout[IDX(xmin,ymax,6)] = fin[IDX(xmin,ymax,8)];
  fout[IDX(xmax,ymax,5)] = fin[IDX(xmax,ymax,7)];
  fout[IDX(xmax,ymin,8)] = fin[IDX(xmax,ymin,6)];

  // for(int x=xmin+1;x<xmax;x++)
  //   {
  //     fout[x][ymin][7] = fin[x][ymin][5];
  //     fout[x][ymin][4] = fin[x][ymin][2];
  //     fout[x][ymin][8] = fin[x][ymin][6];

  //     fout[x][ymax][5] = fin[x][ymax][7];
  //     fout[x][ymax][2] = fin[x][ymax][4];
  //     fout[x][ymax][6] = fin[x][ymax][8];
  //   }
  // for(int y=ymin+1;y<ymax;y++)
  //   {
  //     fout[xmin][y][6] = fin[xmin][y][8];
  //     fout[xmin][y][3] = fin[xmin][y][1];
  //     fout[xmin][y][7] = fin[xmin][y][5];

  //     fout[xmax][y][8] = fin[xmax][y][6];
  //     fout[xmax][y][1] = fin[xmax][y][3];
  //     fout[xmax][y][5] = fin[xmax][y][7];
  //   }
  // fout[xmin][ymin][7] = fin[xmin][ymin][5];
  // fout[xmin][ymax][6] = fin[xmin][ymax][8];
  // fout[xmax][ymax][5] = fin[xmax][ymax][7];
  // fout[xmax][ymin][8] = fin[xmax][ymin][6];
}

      

void computeSquareBounceBack_SURFACE(double ***f, int xmin, int xmax, int ymin, int ymax)
{
  /*Bounce Back on the wet square surface, considering that the square is full of nonfluid matter.
    If you use this function, you should use a streaming & collision procedure that do not apply the streaming and collision to the inside of the square*/
  /*This function is NON OPTMIZED, and its implementation is OBSOLETE*/
  int x0, y0;
  
  /*West side*/
x0 = xmin;
for (int y=ymin+1;y<ymax;y++)
  {
f[x0][y][6] = f[x0+1][y-1][8];
f[x0][y][3] = f[x0+1][y][1];
f[x0][y][7] = f[x0+1][y+1][5];
  }

/*East side*/
x0 = xmax; 
for (int y=ymin+1;y<ymax;y++)
  {
f[x0][y][8] = f[x0-1][y+1][6];
f[x0][y][1] = f[x0-1][y][3];
f[x0][y][5] = f[x0-1][y-1][7];
  }


/*North side*/
y0 = ymax;
for (int x=xmin+1;x<xmax;x++)
  {
f[x][y0][6] = f[x+1][y0-1][8];
f[x][y0][2] = f[x][y0-1][4];
f[x][y0][5] = f[x-1][y0-1][7];
  }

/*South side*/
y0 = ymin;
for (int x=xmin+1;x<xmax;x++)
  {
f[x][y0][8] = f[x-1][y0+1][6];
f[x][y0][4] = f[x][y0+1][2];
f[x][y0][7] = f[x+1][y0+1][5];
  }

/*Top left*/
x0 = xmin; y0 = ymax;
f[x0][y0][6] = f[x0+1][y0-1][8];

/*Top right*/
x0 = xmax; y0 = ymax;
f[x0][y0][5] = f[x0-1][y0-1][7];

/*Bottom right*/
x0 = xmax; y0 = ymin;
f[x0][y0][8] = f[x0-1][y0+1][6];

/*Bottom left*/
x0 = xmin; y0 = ymin;
f[x0][y0][7] = f[x0+1][y0+1][5];

}
	  
	
