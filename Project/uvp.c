#include "uvp.h"
#include "helper.h"
#include <math.h>
#include <stdio.h>

void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V
)
{if(tau>0) 
  {      int i,j;
    double umax=0.0;double vmax=0.0; double mmin;

    for(j=1;j<=jmax;j++)
      for(i=1;i<=imax;i++)
	{umax=fmax(U[i][j],umax);
	  vmax=fmax(V[i][j],vmax);}

    mmin=fmin(dx/fabs(umax),dy/fabs(vmax));
    *dt=tau*fmin(0.5*Re/(1/(dx*dx)+1/(dy*dy)),mmin);
    
  }
}

void calculate_fg(
  double Re,
  double GX,
  double GY,
  double alpha,
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
int **Flag
)
{int i,j;
 
  for(j=1;j<=jmax;j++)
    for(i=1;i<=imax-1;i++)
      { if((Flag[i][j] & 0x10) &&(Flag[i+1][j] & 0x10) !=0)
F[i][j]= U[i][j] + dt*( ((U[i+1][j]-2.0*U[i][j]+U[i-1][j])/(dx*dx) + (U[i][j+1]-2.0*U[i][j]+U[i][j-1])/(dy*dy) )/Re  -1.0/dx*( pow((U[i][j]+U[i+1][j])/2.0,2) - pow((U[i-1][j]+U[i][j])/2.0,2) ) - alpha/dx*( fabs(U[i][j]+U[i+1][j])*(U[i][j]-U[i+1][j])/4.0-fabs(U[i-1][j]+U[i][j])*(U[i-1][j]-U[i][j])/4.0)  - 1.0/dy*((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])/4.0-(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j])/4.0)-alpha/dy*( fabs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])/4.0-(fabs(V[i][j-1]-V[i+1][j-1])*(U[i][j-1]-U[i][j]))/4.0   )+GX);}
  

  for(j=1;j<=jmax-1;j++)
    for(i=1;i<=imax;i++) 
      {if((Flag[i][j] & 0x10) && (Flag[i][j+1] & 0x10) !=0)
G[i][j]= V[i][j] + dt*( ( (V[i+1][j]-2.0*V[i][j]+V[i-1][j])/(dx*dx)+   (V[i][j+1]-2.0*V[i][j]+V[i][j-1])/(dy*dy)  )/Re  - 1.0/dy*(  pow((V[i][j]+V[i][j+1])/2.0,2) - pow((V[i][j-1]+V[i][j])/2.0,2) ) -alpha/dy*(fabs(V[i][j]+V[i][j+1])*(V[i][j]-V[i][j+1])/4.0-fabs(V[i][j-1]+V[i][j])*(V[i][j-1]-V[i][j])/4.0)-1.0/dx*((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])/4.0-(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j])/4.0)-alpha/dx*( fabs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j])/4.0-fabs(U[i-1][j]-U[i-1][j+1])*(V[i-1][j]-V[i][j])/4.0) +GY);}

 
  for(j=1;j<=jmax;j++)
    { F[0][j]=U[0][j];F[imax][j]=U[imax][j];} 
  for(i=1;i<=imax;i++)
    { G[i][0]=V[i][0];G[i][jmax]=V[i][jmax];}

for(j=1;j<=jmax;j++)
    for(i=1;i<=imax;i++){
/*B_N*/
if(Flag[i][j]==1)G[i][j]=V[i][j]; 
/*B_W*/
if(Flag[i][j]==4)F[i-1][j]=U[i-1][j];
/*B_S*/
if(Flag[i][j]==2)G[i][j-1]=V[i][j-1];
/*B_E*/
if(Flag[i][j]==8)F[i][j]=U[i][j];
/*B_NE*/
if(Flag[i][j]==9){F[i][j]=U[i][j];G[i][j]=V[i][j];}
/*B_NW*/
if(Flag[i][j]==5){F[i-1][j]=U[i-1][j];G[i][j]=V[i][j];}
/*B_SE*/
if(Flag[i][j]==10){F[i][j]=U[i][j];G[i][j-1]=V[i][j-1];}
/*B_SW*/
if(Flag[i][j]==6){F[i-1][j]=U[i-1][j];G[i][j-1]=V[i][j-1];}
}
}

void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS
		  )
{int i,j;
 for(j=1;j<=jmax;j++)
    for(i=1;i<=imax;i++)
      RS[i][j]=1.0/dt*((F[i][j]-F[i-1][j])/dx+(G[i][j]-G[i][j-1])/dy  );

}


void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P,int **Flag
)
{int i,j;
for(j=1;j<=jmax;j++)
    for(i=1;i<=imax-1;i++)
if((Flag[i][j] & 0x10) && (Flag[i+1][j] & 0x10) !=0)
      U[i][j]=F[i][j]-dt/dx*(P[i+1][j]-P[i][j]);

 for(j=1;j<=jmax-1;j++)
    for(i=1;i<=imax;i++)
if((Flag[i][j] & 0x10) && (Flag[i][j+1] & 0x10) !=0)
     V[i][j]=G[i][j]-dt/dy*(P[i][j+1]-P[i][j]);

void COMP_TEMP(
double **U,
double **V,
double **TEMP,
int **Flag,
int imax,
int jmax,
double dt,
double dx,
double dy,

double alpha,
double Re,
double Pr)
{
/*initialization and allocation of TEMP in main.c
 TEMP=(double*)malloc(imax*jmax*sizeof(double));*/
double dttxx,dttyy,dutx,dvty;
dttxx=0.0;dttyy=0.0;dutx=0.0;dvty=0.0;

dttxx=(TEMP[i+1][j]-2.0*TEMP[i][j]+TEMP[i-1][j])/(dx*dx);

dttyy=(TEMP[i][j+1]-2.0*TEMP[i][j]+TEMP[i][j-1])/(dy*dy);

dutx=1/dx*( U[i][j]* (TEMP[i][j]+TEMP[i+1][j])/2.0 - U[i-1][j]* (TEMP[i-1][j]+TEMP[i][j])/2.0 )
+gamma/dx*( fabs(U[i][j])* (TEMP[i][j]-TEMP[i+1][j])/2.0  -fabs(U[i-1][j])* (TEMP[i-1][j]-TEMP[i][j])/2.0   );

dvty=1/dy*( V[i][j]* (TEMP[i][j]+TEMP[i][j+1])/2.0 - V[i][j-1]* (TEMP[i][j-1]+TEMP[i][j])/2.0 )
+alpha/dy*( fabs(V[i][j])* (TEMP[i][j]-TEMP[i][j+1])/2.0  -fabs(V[i][j-1])* (TEMP[i][j-1]-TEMP[i][j])/2.0   );

for(j=1;j<=jmax;j++)
for(i=1;i<=imax;i++)
if(Flag[i][j] & 0x10){

TEMP[i][j]=TEMP[i][j]+dt*( 1.0/(Re*Pr)*(dttxx+dttyy)+dvty-dutx);

}

}
}
