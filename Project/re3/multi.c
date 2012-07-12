void multi(  
  double omg,
  double dx,
  double dy,
  int    imax,
  int    jmax,
  double **P,
  double **RS,
  double *res,
  int **Flag,
  char *problem,
  double dp)
{int i,j,numcell;
numcell=0;
itera=0;
 
/* SOR iteration */
  for(i = 1; i <= imax; i++) {
    for(j = 1; j<=jmax; j++) {
      if(Flag[i][j]>=81 && Flag[i][j]<162)
	{P[i][j]
T=v_cycle(T,FR,Nx,Ny);
        }
   }
 }
  /* compute the residual */
  rloc = 0;
  for(i = 1; i <= imax; i++) {
    for(j = 1; j <= jmax; j++) {
      if(Flag[i][j]>=81 && Flag[i][j]<162)
	{
	  rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
	    ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
	}
    }
  }
  rloc = rloc/numcell;
  rloc = sqrt(rloc);
  /* set residual */
  *res = rloc;


ic=(imax-1)/2;Nyc=(Nyf-1)/2;
Dc=Nxc*Nyc;
Tf=GS_multi(Nxf,Nyf,FR,T);
rf=RES_multi(Nxf,Nyf,FR,Tf);
if (Nxf>3)
    rc=restr_multi(Nxc,Nyc,rf);
    ec1=zeros(1,Dc);
    ec2=v_cycle(ec1,rc,Nxc,Nyc);
    ef=interp_multi(Nxf,Nyf,ec2);
    Tf=Tf+ef;
end
Tf=GS_multi(Nxf,Nyf,FR,Tf);


}
