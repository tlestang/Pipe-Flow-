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
	  
	
