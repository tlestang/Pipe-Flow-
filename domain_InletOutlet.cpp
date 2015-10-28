/*Use this function to generate periodic boundary conditions and a pressure gradient beta between inlet (xIn) and outlet (xOut)*/
/*NON OPTIMIZED VERSION*/

void computeDomainInletOutlet(double ***f, int Dx, int Dy, double beta)
{

  int xOut=Dx-1, xIn=0;
  double rhoOut_Avg=0.0, rhoIn_Avg=0.0;
  double rho0 = 1.0;
  double weight, rhoRef;

  /*Compute average densities at Inlet and Outlet*/
  for(int y=1;y<Dy-1;y++)
    {
      for(int k=0;k<9;k++)
	{
	  rhoIn_Avg += f[xIn+1][y][k];
	  rhoOut_Avg += f[xOut-1][y][k];
	}
    }
  rhoOut_Avg /= (Dy-2);
  rhoIn_Avg /= (Dy-2);
  /*Compute reference density*/
  rhoRef = (rhoIn_Avg + rhoOut_Avg)/2.0;
  
  /*Inlet (west wall)*/
  weight = (rho0 + beta*3.0)/rhoOut_Avg;
  for(int y=1;y<Dy-1;y++)
    {
      f[xIn][y][5] = f[xOut-1][y][5]*weight;
      f[xIn][y][1] = f[xOut-1][y][1]*weight;
      f[xIn][y][8] = f[xOut-1][y][8]*weight;
    }
  /*Outlet (east wall)*/
  weight = (rho0 - beta*3.0*(Dx-1))/rhoIn_Avg;			       			       
  for(int y=1;y<Dy-1;y++)
    {
      f[xOut][y][6] = f[xIn+1][y][6]*weight;
      f[xOut][y][3] = f[xIn+1][y][3]*weight;
      f[xOut][y][7] = f[xIn+1][y][7]*weight;
    }
}
