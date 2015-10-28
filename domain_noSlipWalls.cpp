void computeDomainNoSlipWalls_BB(double ***fout, double ***fin, int Dx, int Dy)
{
  int y0;
  for(int x=0;x<Dx;x++)
    {

      /*North boundary*/
      fout[x][Dy-1][4] = fin[x][Dy-1][2];
      fout[x][Dy-1][8] = fin[x][Dy-1][6];
      fout[x][Dy-1][7] = fin[x][Dy-1][5];

      /*South boundary*/
      fout[x][0][2] = fin[x][0][4];
      fout[x][0][5] = fin[x][0][7];
      fout[x][0][6] = fin[x][0][8];
    }
}
