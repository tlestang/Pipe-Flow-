void initializePopulations(double ***fin, int Dx, int Dy)
{
  double w[9]={4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
  for(int x=0;x<Dx;x++)
    {
      for(int y=0;y<Dy;y++)
	{
	  for(int k=0;k<9;k++)
	    {
	      fin[x][y][k] = w[k];
	    }
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
