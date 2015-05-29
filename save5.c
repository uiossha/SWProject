/*=========================================
  This program sets up the matrix equation [A]{f}={Q} which
  results from finite volume discretization of the
  Poissson equation in 2D using Cartesian grid and
  central difference approximation of the second derivative.
  Boundary conditions are of Dirichlet type on east and north
  boundary: f(1,y) = 0, f(x,1) = 0., and of Neumann type on
  west and south boundary (zero gradient normal to boundary).
  The equation can be solved by a variety of iterative
  C     solvers.
  ===============================================*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "data.h"
#define min(a,b)\
  ({__typeof__ (a) _a = (a);\
    __typeof__ (b) _b = (b);\
    _a < _b ? _a : _b; })

struct space{
  int ni;     /*number of points on line i*/
  int nj;     /*number of points on line j*/
  int njm;
  int nim;    /*number of meshs on line i*/
  double ae[NXY];/*area's east*/
  double as[NXY];/*area's south*/
  double aw[NXY];/*area's west*/
  double an[NXY];/*area's north*/
  double ap[NXY];/*area's points*/
  double fi[NXY];/*Field*/
  double q[NXY]; 
  double xc[NX];/*x cells*/
  double yc[NY];/*y cells*/
  int li[NX]; /*list*/
  int x[NX];
  int y[NY];
};

struct vars{   /*variables*/
  int maxit;  /*maximum iteration*/
  double resmax;/*max of residue*/   
  double errnor;/*error norm*/ 
  double res0;  /*residueal sum*/
};

typedef struct space myspace;
typedef struct vars myvars;

void GSS(myspace s1, myvars v1, double alfa);

void LSOL(myspace s1, myvars v1);

void ADI(double beta,myspace s1, myvars v1);

int main(void){
  
  myspace s1;
  myvars v1;
  int n = 25;
  int m = 25;
  double dx, dy;  /* dimension x and y*/
  
  double dim[NXY][NX][NY];
  
  /*Maximum number of allowed iteration and normalized residuel sum*/
  printf("Enter: max iteration\n");
  scanf("%d", &v1.maxit);
  printf("Your max iteration is %d\n", v1.maxit);
  printf("Enter: max of residual\n");
  scanf("%lf", &v1.resmax);
  printf("Your max of residual is %lf\n", v1.resmax);
  /*DEFINE GRID IN X-DIRECTION*/
  
  
  double xmax, exx, nicv;
  printf("Enter: max of x\n");
  scanf("%lf", &xmax);
  printf("Your xmax is %lf\n", xmax);
  printf("Enter: exx\n");
  scanf("%lf", &exx);
  printf("Your exx is %lf\n", exx);
  printf("Enter: nicv\n");
  scanf("%lf", &nicv);
  printf("Your nicv is %lf\n", nicv);
  
  s1.ni = n+2;
  s1.nim = s1.ni-1;
  if(exx == 1)
    {dx=xmax/n;}
  s1.x[1]=0;
  for(int i = 2; i<=s1.nim; i++){
    s1.x[i] = s1.x[i-1]+dx;
    dx=dx*exx;
  }
  
  /*DEFINE GRID IN Y-DIRECTION*/
  
  double ymax, exy, njcv;
  printf("Enter: ymax\n");
  scanf("%lf", &ymax);
  printf("Your ymax is %lf\n", ymax);
  printf("Enter: exy\n");
  scanf("%lf", &exy);
  printf("your exy is %lf\n", exy);
  printf("Enter: njcv\n");
  scanf("%lf", &njcv);
  printf("Your njcv is %lf\n", njcv);
  
  s1.nj=m+2;
  s1.njm=s1.nj-1;
  if(exy == 1)
    {dy=ymax/m;}
  else{
    dy=ymax*(1.0-exy)/(1.0-pow(exy, m));
  }
  s1.y[1]=0.0;
  for(int j = 2;j<=s1.njm;j++){
    s1.y[j]=s1.y[j-1]+dy;
    dy=dy*exy;
  }
  /*COORDINATES OF CELL CENTERS*/
  
  s1.xc[1]=s1.x[1];
  s1.xc[s1.ni]=s1.x[s1.nim];
  for(int i = 2;i<=s1.nim;i++)
    {s1.xc[i]=0.5*(s1.x[i]+s1.x[i-1]);}
  
  s1.yc[1]=s1.y[1];
  s1.yc[s1.nj]=s1.y[s1.njm];
  for(int j = 2; j<=s1.njm;j++){
    s1.yc[j]=0.5*(s1.y[j]+s1.y[j-1]);}
  
  /*WORKING ARRAY FOR CONVERTING 2D INDICES TO 1D*/
  
  int nij = s1.ni*s1.nj;
  for(int i = 1;i<=s1.ni;i++){
    s1.li[i] = (i-1)*s1.nj;
  }
  
  /*INITIALIZE FIELD VALUES(ZERO)*/
  
  for(int ij=1;ij<=s1.nj;ij++){
    s1.fi[ij]=0;
  }
  
  /*From the earlier version: Laplace equation with exact solution
    
    C     Fi = x * y when Dirichlet boundary condition is applied at all
    
    C     boundaries (linear profiles, CDS gives exact solution on any grid).
    
    C     In that case, FI=0 at south and west boundary remains after 
    
    C     initialization, and the boundary values at east and north boundaries
    
    C     must be prescribed:*/
  
  
  for(int i=2; i<= s1.ni; i++){
    s1.fi[s1.li[i]+s1.nj] = s1.xc[i]*s1.yc[s1.nj];}
  for(int j =2 ; j<= s1.njm; j++){
    s1.fi[s1.li[s1.ni]+j]=s1.yc[j]*s1.xc[s1.ni];}
  
  /*C
    C.....CALCULATE ELEMENTS OF MATRIX [A] (For Laplace equation, Q(IJ)=0.;
    
    C     for Poisson equation, the source term should be chosen so that
    
    C     the sum of all sources = 0 if Neumann conditions are applied on 
    
    C     all boundaries.). Here, the sourse term is a sine function...
    
    C*/
  
  
  for(int i = 2; i<=s1.nim; i++){
    for(int j=2; j<=s1.njm; j++){
      int ij = s1.li[i]+j;
      s1.ae[ij]=-(s1.y[j]-s1.y[j-1])/(s1.xc[i+1]-s1.xc[i]);
      s1.an[ij]=-(s1.x[i]-s1.x[i-1])/(s1.yc[j+1]-s1.xc[j]);
      s1.aw[ij]=-(s1.y[j]-s1.y[j-1])/(s1.xc[i]-s1.xc[i-1]);
      s1.as[ij]=-(s1.x[i]-s1.x[i-1])/(s1.yc[j]-s1.yc[j-1]);
      s1.ap[ij]=-(s1.ae[ij]+s1.aw[ij]+s1.an[ij]+s1.as[ij]);
      double vol=(s1.x[i]-s1.x[i-1]*(s1.y[j]-s1.y[j-1]));
      s1.q[ij]=sin(2*M_PI*s1.xc[i]/xmax)*sin(2*M_PI*s1.yc[j]/ymax)*vol;
    }
  }
  /*C
    C.....NEUMANN CONDITION AT BOUNDARIES X=0 & Y=0: quadratic extrapolation
    C     from inside is used, boundary value is expressed through a weighted
    C     average of two inner nodes, and the product of the coefficient 
    C     multiplying boundary noal value leads then to the modification of
    C     the coefficients for two inner nodes involved. The Neumann b.c. is
    C     thus implicitly taken into account; the boundary values have to be
    C     computed in each iteration if multigrid method is used, otherwise
    C     only at the end.
    C*/
  
  int i=2;
  double quo = pow((s1.xc[3]-s1.xc[1]),2)-pow((s1.xc[2]-s1.xc[1]),2);
  for(int j=2;j<=s1.njm;j++){
    int ij=s1.li[i]+j;
    s1.ap[ij]=s1.ap[ij]+s1.aw[ij]*pow((s1.xc[3]-s1.xc[1]),2)/quo;
    s1.ae[ij]=s1.ae[ij]-s1.aw[ij]*pow((s1.xc[2]-s1.xc[1]),2)/quo;
    s1.aw[ij]=0.0;
  }
  int j=2;
  quo=pow((s1.yc[3]-s1.yc[1]),2)-pow((s1.yc[2]-s1.yc[1]),2);
  for(int i = 2;i<=s1.nim;i++){
    int ij =s1.li[i]+j;
    s1.ap[ij]=s1.ap[ij]+s1.as[ij]*pow((s1.yc[3]-s1.yc[1]),2)/quo;
    s1.an[ij]=s1.an[ij]-s1.as[ij]*pow((s1.yc[2]-s1.yc[1]),2)/quo;
    s1.as[ij]=0;
  }
  /*C
    C.....CHOOSE THE SOLVER
    C*/
  
  printf("Enter solver ID:\n");
  printf("1 - Gauss-Seidel solver\n");
  printf("2 - Line-by-Line TDMA\n");
  printf("3 - Stone SIP solver\n");
  printf("4 - Conjugate gradient solver\n");
  printf("5 - ADI\n");
  printf("6 - Multigrid solver with GS\n");
  printf("7 - Multigrid solver with SIP\n");
  printf("8 - Multigrid solver with ICCG\n");
  printf("  \n");
  int isol;
  scanf("%d", &isol);	
  double alfa, beta;
  if(isol == 3 || isol == 5 || isol == 7){
    printf("Enter: alfa, beta (relevant only for SIP or ADI\n");
    printf("        \n");
    scanf("%lf%lf", &alfa, &beta);
    printf("Your alfa and beta is %lf , %lf\n", alfa, beta);
  }
  
  /*
    C
    C.....CALCULATE INITIAL ERROR NORM (only for the laplace equation, where
    C     the solution is known, Fi = x*y)
    C
    c      ERRNOR=0.
    c      DO I=2,NIM
    c        DO J=2,NJM
    c          IJ=LI(I)+J
    c          ERRNOR=ERRNOR+ABS(FI(IJ)-XC(I)*YC(J))
    c        END DO
    c      END DO
    c      WRITE(8,*) '  Initial Error norm:  ERR = ',ERRNOR
    C
    C.....SOLVE EQUATION SYSTEM
    C*/
  switch(isol){
  case 1:
    printf("You have chosen the ");
    printf("Gauss-Seidel solver\n");
    GSS(s1, v1, alfa);
    break;
  case 2:
    printf("You have chosen the ");
    printf("Line-by-Line TDMA\n");
    /* LSOL(s1, v1, alfa);*/
    break;
  case 3:
    printf("You have chosen the ");
    printf("Stone SIP solver\n");
    /*CALL SIPSOL(FI)*/
    break;
  case 4:
    printf("You have chosen the ");
    printf("Conjugate gradient solver\n");
    /* CALL CGS(FI)*/
    break;
  case 5:
    printf("You have chosen the ");
    printf("ADI\n");
    /*CALL ADI(BETA,FI)*/
    break;
  case 6:
    printf("You have chosen the ");
    printf("Multigrid solver with GS\n");
    /* CALL MG(NI,NJ,ISOL,FI)*/
    break;
  case 7:
    printf("You have chosen the ");
    printf("Multigrid solver with SIP\n");
    /* CALL MG(NI,NJ,ISOL,FI)*/
    break;
  case 8:
    printf("You have chosen the ");
    printf("Multigrid solver with ICCG\n");
    /* CALL MG(NI,NJ,ISOL,FI)*/
    break;
  default:
    break;
    
  }
  /*
    C
    C.....PRINT ERROR NORM AND SOLUTION FIELD (applies only to Laplace
    C     equation with Dirichlet b.c.; for the Poisson equation, the exact
    C     solution is not known...
    C*/
  double err = 0;
  for(int i = 2; i<= s1.nim; i++){
    for(int j = 2; j<= s1.njm; j++){
      int ij = s1.li[i]+j;
      err = err + fabs(s1.fi[ij]-s1.xc[i]*s1.yc[j]);
    }
  }
  err = err/v1.errnor;
  printf("Final error norm: err = %f",err);
  /*
    IF(ISOL.EQ.1) WRITE(8,*) '  GAUSS-SEIDEL SOLVER'
    IF(ISOL.EQ.2) WRITE(8,*) '  LINE-BY-LINE (ADI) TDMA SOLVER'
    IF(ISOL.EQ.3) WRITE(8,*) '  SIP SOLVER'
    IF(ISOL.EQ.4) WRITE(8,*) '  ICCG SOLVER'
    IF(ISOL.EQ.5) WRITE(8,*) '  ADI SOLVER'
    IF(ISOL.EQ.6) WRITE(8,*) '  MG-GS SOLVER'
    IF(ISOL.EQ.7) WRITE(8,*) '  MG-SIP SOLVER'
    IF(ISOL.EQ.8) WRITE(8,*) '  MG-ICCG SOLVER'
  */
   printf("ni = %d, nj = %d, 0 , fi = %fl", s1.ni,s1.nj,s1.fi[0]);
}

/*C##############################################################
  SUBROUTINE GSS(FI)
  C##############################################################
  C     This routine contains the Gauss-Seidel solver
  C==============================================================*/
void GSS(myspace s1, myvars v1, double alfa){
  for(int n = 1; n<=v1.maxit;n++){
    /*	C.....CALCULATE "FALSE" RESIDUAL AND UPDATE VARIABLE*/
    for(int i = 2; i<=s1.nim; i++){
      for(int ij = s1.li[i] + 2 ; ij<= s1.li[i] + s1.njm; ij++){
	double res = s1.q[ij]-s1.ap[ij]*s1.fi[ij]-	       \
	  s1.ae[ij]*s1.fi[ij+s1.nj]-s1.aw[ij]*s1.fi[ij-s1.nj]- \
	  s1.as[ij]*s1.fi[ij-1]-s1.an[ij]*s1.fi[ij+1];
	s1.fi[ij]=s1.fi[ij]+res/s1.ap[ij];
      }
    }
    
    /*C.....CHECK CONVERGENCE*/
    
    double resn=0.0;
    for(int i = 2; i<=s1.nim; i++){
      for(int j = 2; j<=s1.njm; j++){
	int ij = s1.li[i]+j;
	double res = s1.q[ij]-s1.ap[ij]*s1.fi[ij]- \
	  s1.ae[ij]*s1.fi[ij+s1.nj]-s1.aw[ij]*s1.fi[ij-s1.nj]-\
	  s1.as[ij]*s1.fi[ij-1]-s1.an[ij]*s1.fi[ij+1];
	resn=resn+fabs(res);
      }
    }
    
    if(n == 1){
      v1.res0=resn;
    }
    double rsm = resn/v1.res0;
    printf("Sweep, rsm = %lf \n", rsm);
    if(rsm<v1.resmax){
      return;
    }
  }
  return ;
}
