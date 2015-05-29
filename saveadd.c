/*

C###########################################################     
      SUBROUTINE LSOL(FI)
C###########################################################
C     In this solver the TDMA algorithm is applied line-by-line
C     alternately along J and I lines.
C     Warning: arrays A() and C() should be dimensioned
C     as the maximum of NX and NY!
C===========================================================
*/

void LSOL(int NX, int NY, int NI, int NJ, int NIM, int NJM, int MAXIT, int LI[], double RESMAX, double ALFA, double ERRNOR, double RES0, double X[], double Y[], double AE[], double AW[], double AN[], double AS[], double AP[], double Q[],double FI[], double XC[], double YC[]){

      double A[NX],C[NX];
      /*DATA A,C /NX*0.,NX*0./ */
      for(int N = 1; N <= MAXIT; N++){    
/*C.....SOLVE WITH "TDMA" ALONG  J-LINES*/
		for(int i = 2; i<= NIM; i++){
			for(int j = 2; j<= NJM; j++){
				int ij = LI[i]+j;
            	double APR = (1.0)/(AP[ij]-AS[ij]*A[j-1]);
            	A[j]=AN[ij]*APR;
            	C[j]=(Q[ij]-AE[ij]*FI[ij+NJ]-AW[ij]*FI[ij-NJ]-AS[ij]*C[j-1])*APR;
            }
            
        	for(int j =NJM; j>=2; j--){
        		int ij = LI[i]+j;
        		FI[ij] = C[j]-A[j]*FI[ij+1];
        	}
        }
        
/*C.....SOLVE WITH "TDMA" ALONG  I-LINES*/
		for(int j = 2; j<=NJM; j++){
			for(int i = 2; i<=NIM; i++){
				int ij = LI[i]+j;
				double APR= (1.0)/(AP[ij]-AW[ij]*A[i-1]);
				A[i] = AE[ij]*APR;
				C[i]=(Q[ij]-AN[ij]*FI[ij+1]-AS[ij]*FI[ij-1]-AW[ij]*C[i-1])*APR;
			}
			
			for(int i=NIM; i>=2; i--){
				int ij = LI[i]+j;
				FI[ij]=C[i]-A[i]*FI[ij+NJ];
			}
		}
/*.....CHECK CONVERGENCE OF INNER ITERATIONS*/

        double RESN=0.0;
        for(int i = 2; i<=NIM; i++){
			for(int j = 2; j<=NJM; j++){
				int ij = LI[i]+j;
				double RES = Q[ij]-AP[ij]*FI[ij]-AE[ij]*FI[ij+NJ]-AW[ij]*FI[ij-NJ]-AS[ij]*FI[ij-1]-AN[ij]*FI[ij+1];
				RESN = RESN + abs(RES);
			}
		}
		    
		if(N == 1){RES0=RESN;}
		double RSM = RESN/RES0;
        printf("SWEEP, RSM = %lf \n", RSM);
        if(RSM<RESMAX){
        	return;
        }
    }
    return;
}

/*
C###########################################################     
      SUBROUTINE ADI(BETA,FI)
C###########################################################
C     This is the standard ADI solver, as described in Sect.
C     5.3.5. The original version of this code was written
C     by Joel H. Ferziger, Stanford University, 1995.
C===========================================================
*/

void ADI(double BETA, int NX, int NY, int NI, int NJ, int NIM, int NJM, int MAXIT, int LI[], double RESMAX, double ALFA, double ERRNOR, double RES0, double X[], double Y[], double AE[], double AW[], double AN[], double AS[], double AP[], double Q[],double FI[], double XC[], double YC[]){
	  int NXY = NX*NY;
      double AI[NX], CI[NX], AJ[NY], CJ[NY], API[NX], QI[NX], APJ[NY], QJ[NY], FF[NXY]; 

/*C.....CALCULATE ADI COEFFICIENTS*/
	
	for(int i =1; i<= NI; i++){
		for(int j = 1; j<=NJ; j++){
			int ij = LI[i]+j;
		}
	}
	/*.....START ITERATION LOOP*/
	for(int N = 1; N<=MAXIT; N++){
		
/*C.....SOLVE WITH "TDMA" ALONG  I-LINES*/
		for(int j = 2; j<=NJM; j++){
			for(int i = 2; i<=NIM; i++){
				int ij = LI[i]+j;
				API[i]=BETA-AE[ij]-AW[ij];
				QI[i]=BETA+AN[ij]+AS[ij];
			}
			for (int i = 2; i<=NIM; i++){
			int ij = LI[i]+j;
			double APR = (1.0)/(API[i]-AW[ij]*AI[i-1]);
			AI[i]=AE[ij]*APR;
			CI[i]=(Q[ij]+QI[i]*FI[ij]-AN[ij]*FI[ij+1]-AS[ij]*FI[ij-1]-AW[ij]*CI[i-1])*APR;
		}
		for(int i = NIM; NIM>=2; i--){
			int ij = LI[i]+1;
			FF[ij] = CI[i] - AI[i]*FF[ij+NJ];
		}
	}
	    
/*C.....SOLVE WITH "TDMA" ALONG  J-LINES*/
	for(int i = 2; i<=NIM; i++){
		for(int j = 2; j<=NJM; j++){
			int ij = LI[i] +j ;
			APJ[j]=BETA-AN[ij]-AS[ij];
			QJ[j]=BETA+AE[ij]+AW[ij];
		}
		for(int j = 2; j<=NJM; j++){
			int ij = LI[i]+j;
			double APR = (1.0)/(APJ[j]-AS[ij]*AJ[j-1]);
			AJ[j]=AN[ij]*APR;
			CJ[j]=(Q[ij]+QJ[j]*FF[ij]-AE[ij]*FF[ij+NJ]-AW[ij]*FF[ij-NJ]-AS[ij]*CJ[j-1])*APR;
		}
		for (int j=NJM; j>=2; j--){
			int ij = LI[i]+j;
			FI[ij] = CJ[j]-AJ[j]*FI[ij+1];
		}
	}
	
/*C.....CHECK CONVERGENCE OF INNER ITERATIONS*/

    	double RESN=0.0;
    	for(int i = 2; i<= NIM; i++){
    		for(int j = 2;j<=NJM; j++){
    			int ij = LI[i]+j;
    			double RES = Q[ij]-AP[ij]*FI[ij]-AE[ij]*FI[ij+NJ]-AW[ij]*FI[ij-NJ]-AS[ij]*FI[ij-1]-AN[ij]*FI[ij+1];
    			RESN = RESN + abs(RES);
    		}
    	}
    	
    	if(N == 1){RES0=RESN;}
    	double RSM = RESN/RES0;
    	printf("SWEEP, RSM = %lf \n", RSM);
        if(RSM<RESMAX){
        	return;
        }
    }
    return;
}
    	/*

C
C#############################################################
      SUBROUTINE SIPSOL(FI)
C#############################################################
C     This is the ILU solver after Stone; see Sect. 5.3.4 for
C     a description of the algorithm.
C
C     M. Peric, Institut fuer Schiffbau, Hamburg, 1995
C=============================================================
*/
void SIPSOL(double BETA, int NX, int NY, int NI, int NJ, int NIM, int NJM, int MAXIT, int LI[], double RESMAX, double ALFA, double ERRNOR, double RES0, double X[], double Y[], double AE[], double AW[], double AN[], double AS[], double AP[], double Q[],double FI[], double XC[], double YC[]){

	  int NXY=NX*NY;
      double LW,LS,LPR;
      double UN[NXY],UE[NXY],RES[NXY], \
               LW2[NXY],LS2[NXY],LPR2[NXY];

/*C.....CALCULATE ELEMENTS OF [L] AND [U] MATRICES*/

	for(int i =2; i<=NIM; i++){
		for(int ij = LI[i]+2; ij<=LI[i]+NJM; ij++){
			LW2[ij]=AW[ij]/((1.0)+ALFA*UN[ij-NJ]);
			LS2[ij]=AS[ij]/((1.0)+ALFA*UE[ij-1]);
			double P1=ALFA*LW2[ij]*UN[ij-NJ];
			double P2=ALFA*LS2[ij]*UE[ij-1];
			LPR2[ij]=(1.0)/(AP[ij]+P1+P2-LW2[ij]*UE[ij-NJ]-LS2[ij]*UN[ij-1]);
			UN[ij]=(AN[ij]-P1)*LPR2[ij];
			UE[ij]=(AE[ij]-P2)*LPR2[ij];
		}
	}
	
    
/*C.....CALCULATE RESIDUAL AND AUXILLIARY VECTORS; INNER ITERATION LOOP*/

	for(int N=1; N=MAXIT; N++){
		
		int RESN=0.0;
		
		for	(int i = 2; i<=NIM; i++){
			for(int ij=LI[i]+2; ij<=LI[i]+NJM; ij++){
				RES[ij]=Q[ij]-AP[ij]*FI[ij]-AN[ij]*FI[ij+1]- \
					AS[ij]*FI[ij-1]-AN[ij]*FI[ij+NJ]-AW[ij]*FI[ij-NJ];
					RESN=RESN+abs(RES[ij]);
					RES[ij]=(RES[ij]-LS2[ij]*RES[ij-1]-LW2[ij]*RES[ij-NJ])*LPR2[ij];
				}
			}
			if(N == 1){RES0 = RESN;}
		
/*C.....CALCULATE INCREMENT AND CORRECT VARIABLE*/
		for(int i = NIM; i>=2; i--){
			for(int ij=LI[i]+NJM; ij<=LI[i]+2; ij--){
				RES[ij]=RES[ij]-UN[ij]*RES[ij+1]-UE[ij]*RES[ij+NJ];
				FI[ij]=FI[ij]+RES[ij];
			}
		}
		

/*C.....CONVERGENCE CHECK*/

		double RSM = RESN/(RES0+1.e-20);
		printf("SWEEP, RSM = %lf \n", RSM);
        if(RSM<RESMAX){
        	return;
        }
    }
    return;
}

/*
C##############################################################
      SUBROUTINE CGS(FI)
C##############################################################
C     This is a pre-conditioned conjugate gradient solver for
C     symmetric matrices (e.g., pressure or pressure-correction
C     equation, heat conduction, etc.). A 3D version is included
C     in LAPL3D.F file. The original code was written by Ismet
C     Demirdzic, Masinski Fakultet, Sarajevo, 1987 (see Sect.
C     5.3.6 for a description of the algorithm).
C===============================================================
*/ 
void CGS(double BETA, int NX, int NY, int NI, int NJ, int NIM, int NJM, int MAXIT, int LI[], double RESMAX, double ALFA, double ERRNOR, double RES0, double X[], double Y[], double AE[], double AW[], double AN[], double AS[], double AP[], double Q[],double FI[], double XC[], double YC[]){

	int NXY=NX*NY;
	double PK[NXY], ZK[NXY], D[NXY], RES[NXY];
		
/*C.....CALCULATE INITIAL RESIDUAL VECTOR*/

	RES0=0;
	for(int i = 2; i<=NIM; i++){
		for(int ij = LI[i]+2; ij <= LI[i]+NJM; ij++){
			RES[ij]=Q[ij]-AP[ij]*FI[ij]-AE[ij]*FI[ij+NJ]- \
			AW[ij]*FI[ij-NJ]-AN[ij]*FI[ij+1]-AS[ij]*FI[ij-1];
			RES0=RES0+abs(RES[ij]);
		}
	}
	

/*C.....PRECONDITIONING MATRIX DIAGONAL*/

	for(int i = 2; i<=NIM; i++){
		for(int ij=LI[i]+2; ij<=LI[i]+NJM; ij++){
			D[ij]=(1.0)/(AP[ij]-D[ij-NJ]*pow(AW[ij],2)-D[ij-1]*pow(AS[ij], 2));
		}
	}
	
/*C.....CALCULATION OF  ZK; INNER ITERATION LOOP*/

	double S0 = 1.e20;
	for(int n = 1; n<=MAXIT; n++){
		
/*C.....FORWARD SUBSTITUTION*/
		for(int i = 2; i<=NIM; i++){
			for(int ij=LI[i]+2; ij<=LI[i]+NJM; ij++){
				ZK[ij]=(RES[ij]-AW[ij]*ZK[ij-NJ]-AS[ij]*ZK[ij-1])*D[ij];
			}
		}
	

		for (int i = 2; i<=NIM; i++){
			for(int ij= LI[i]+2; ij<=LI[i]+NJM; ij++){
				ZK[ij]=ZK[ij]/(D[ij]+1.e-20);
			}
		}
		
/*C.....BACKWARD SUBSTITUTION*/

	double SK = 0.0;
	for (int i = NIM; i>=2; i--){
		for(int ij = LI[i]+NJM; ij>=LI[i]+2; ij--){
			ZK[ij]=(ZK[ij]-AE[ij]*ZK[ij+NJ]-AN[ij]*ZK[ij+1])*D[ij];
			SK = SK + RES[ij]*ZK[ij];
		}
	}
	
/*C.....CALCULATE BETA AND NEW SEARCH VECTOR*/

	double BET = SK/S0;
	for(int i = 2; i<=NIM; i++){
			for(int ij= LI[i]+2; ij<=LI[i]+NJM; ij++){
				PK[ij] = ZK[ij]+BET*PK[ij];
			}
		}
		
		
/*C.....CALCULATE SCALAR PRODUCT (PK . A PK) AND ALPHA (A PK OVERWRITES ZK)*/

		double S2 = 0.0;
		for(int i = 2; i<=NIM; i++){
			for(int ij = LI[i]+2; ij<=LI[i]+NJM; ij++){
				ZK[ij]= AP[ij]*PK[ij]+AE[ij]*PK[ij+NJ]+AW[ij]*PK[ij-NJ]+\
				AN[ij]*PK[ij+1]+AS[ij]*PK[ij-1];
				S2 = S2+PK[ij]*ZK[ij];
			}
		}
	
		double ALF= SK/S2;
		
		
/*C.....CALCULATE NEW RESIDUAL AND UPDATE VARIABLE*/
		double RESN = 0;
		for(int i = 2; i<=NIM; i++){
			for(int ij = LI[i]+2; ij<= LI[i]+NJM; ij++){
				FI[ij]=FI[ij]+ALF*PK[ij];
				RES[ij]=RES[ij]-ALF*ZK[ij];
				RESN=RESN+abs(RES[ij]);
			}
		}
		S0=SK;
		
/*C....CHECK CONVERGENCE*/
	double RSM = RESN/(RES0+1.e-20);
	printf("SWEEP, RSM = %lf \n", RSM);
        if(RSM<RESMAX){
        	return;
        }
    }
    return;
}
	
	void PRINT2(int NI, double NJ, double IJST, double FI[], char NAME, int NX){
		/*
C#########################################################
      SUBROUTINE PRINT(NI,NJ,IJST,FI,NAME)
C#########################################################
      PARAMETER (NX=162,NY=162,NXY=NX*NY,NXYA=NXY+NXY/2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FI(NXYA),LI(NX)
      CHARACTER*6 NAME
*/
	
	int LI[NX];

	for(int i = 1; i<=NI; i++){
		LI[i]=IJST+(i-1)*NJ;
	}
	
	printf("NAME\n");
	int IE=1;
	double NL = NI/12 +1;
	if(NI%12 == 0){NL = NI/12;}
	for(int L = 1; L<= NL ;L++){
		int IS = IE;
		IE = min(NI, IS+11);
		int i = IS;
		printf("i= %d , IS = %d, IE  = %d", i, IS, IE);
		for(int j = NJ; j>=1; j--){
			printf("%lf, %lf, %lf\n", FI[LI[i]+j], i=IS, IE);
		}
	}
	return;
}
	


