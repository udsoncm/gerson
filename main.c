#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

int periodic(int N, double *diag, double *evals, double *evecs) {
  int i;
  for (i=0; i<N*N; i++)
    evecs[i] = 0.0;
    evals[i] = 0.0;
  for (i=0; i<N; i++)
    evecs[i*N+i] = diag[i];
  for (i=0; i<N-1; i++)
    evecs[i*N+i+1] = -1.0;
  evecs[0*N+N-1] = -1.0;

  LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'L', N, evecs, N, evals);

  return(0);
}

int hardwall(int N, double *diag, double *subd, double *evals, double *evecs, int *ifail) {
  int Nf;
  int Nemin = 1, Nemax = N;
  double Emin = 0.0, Emax = 0.0;
  double abstol = 1E-4;

  LAPACKE_dstevx(LAPACK_COL_MAJOR, 'V', 'A', N, diag, subd, Emin, Emax, Nemin, Nemax, abstol, &Nf, evals, evecs, N, ifail);

  return(0);
}

/* int density(int N, double *L, double *dx, double *Efermi, double *evals, double *evecs, double *dens){ */
/*   int i; */
/*   int j; */
/*   double x; */
/*   FILE *fp1 = fopen("density.dat", "w"); */
    
/*     for ( i = 0; i < N; i++) { */
/*         for (j=0; j<N; j++) */
/*             dens[j] = 0.0; */
/*             while (evals[i] < Efermi) { */
/*             dens[j] = dens[j] + pow(evecs[i*N+j], 2.0); */
/* 	    } // while */
/*     } */
/*     for (j=0; j<N; j++) { */
/*       x = -L/2.0 + j*dx; */
/*       fprintf(fp1, "%g %g\n", x, dens[j]); */
/*     } // for j */
/*   return(0); */
/* } */

int main() {
  int i;
  int N = 1000; // number of points
  double dx = 1.0; // nm
  // x = -L/2 + i*dx
  double L = (N-1)*dx, x; // nm
  double wx = L/10.0; // nm
  double meff = 0.067; // GaAs
  double ep0 = 12.35; // GaAs Dielectric constant
  double Ry = 13.6056917253e3; // meV
  double a0 = 0.0529177208319; // nm
  double h2m = 2.0*Ry*a0*a0/meff/2.0; // h^2/2m [meV nm^2] => Ry = hbar^2/2/m/a0^2
  double Eunit = h2m/dx/dx;
  double Ry_ef = ( Ry * meff ) / (ep0 * ep0);
  double a0_ef = (ep0/meff)*a0;  
  double V0 = 50.0/Eunit; // meV
  double lambda = 1.0;  // hartree potential strenght 
  double *v = (double *) calloc(N, sizeof(double));
  double *dens = (double *) calloc(N, sizeof(double)); // electron density
  double *poth = (double *) calloc(N, sizeof(double)); // electron density
  double *potsc = (double *) calloc(N, sizeof(double)); // electron density
  double *evalsold = (double *) calloc(N, sizeof(double)); // electron density
  FILE *fp = fopen("ldos.dat", "w");
  FILE *fp1 = fopen("density.dat", "w");
  FILE *fp2 = fopen("density.dat", "w");
  double Efermi = 50.0/Eunit; // Fermi level
  double alpha = 0.07;
  double densx;
  double deltaE;
  double deltaE1;
  int j;
  int k;
  int nit = 1000; // number of self-consistent cycles. 

  

  // LAPACK
  double *diag = (double *) calloc(N, sizeof(double));
  double *subd = (double *) calloc(N-1, sizeof(double));
  double *evals = (double *) calloc(N, sizeof(double));
  double *evecs = (double *) calloc(N*N, sizeof(double));
  int *ifail = (int *) calloc(N, sizeof(int));
  double *Vh = (double *) calloc(N, sizeof(double));
  // ==============

  for (i=0; i<N; i++) {
    x = -L/2.0 + i*dx;
    v[i] = V0*exp(-pow(x/wx, 2.0)/2.0); // gaussian barrier
    diag[i] = 2.0 + v[i];
  }
  for (i=0; i<N-1; i++){
    subd[i] = -1.0;
  }

  for (i=0; i<N; i++){
    poth[i] = 0.0;
    potsc[i] = 0.0;
    dens[i] = 0.0;
    evalsold[i] = 1.0/Eunit;
  }


  for (k = 0; k < 1000; k++){ // begin self-consistent cycle 

    for ( i = 0; i < N; i++) { // generating sc potential 
      poth[i] = (2.0 * lambda * dens[i] * Ry_ef)/Eunit;
      //     potsc[i] = v[i];
      potsc[i] = alpha * (v[i] + poth[i]) + (1.0-alpha) * potsc[i];
      diag[i] = 0.0;
      evals[i] = 0.0;
      diag[i] = 2.0 + potsc[i];
    }

   // diagonalizing the Hamiltonian

  //  hardwall(N, diag, subd, evals, evecs, ifail);
    periodic(N, diag, evals, evecs);

    deltaE = evals[1]-evalsold[1];
    deltaE1 = fabs(deltaE);

    if (deltaE1 < 1.0E-8 ){
      printf("program converged   " );
        break; 
    }
    else {
         densx = 0.0;    
	 for (i=0; i< N; i++){
	  evalsold[i] = evals[i]; 
          densx = densx + dens[i];
        
	  }
	 printf("%d %g %g\n", k, densx, evals[1] * Eunit);
    }

    for (i=0; i<N; i++){
        dens[i] = 0.0;
    }

   // calculating electron density
      for ( i = 0; i < N; i++) {
       if (evals[i] <= Efermi) {
	  for (j=0; j<N; j++) { 
              dens[j] += 2.0 * pow(evecs[i*N+j], 2.0);
	  }		
       }   // while
   }



  } // end SC cycle

     for (j=0; j<N; j++) {
       x = -L/2.0 + j*dx;
       fprintf(fp1, "%g %g\n", x, dx*dens[j]);
    } // for j


  //for (i=0; i<N; i++) fprintf(fp, "%g\n", evals[i]*Eunit);
  // LDOS(E) = sum_{i~E} |psi(x)|^2
  // file cols: x E LDOS
  int step, steps = 1000;
  double Emin = evals[0], Emax = evals[N-1];
  double Ei, dE = (Emax-Emin)/steps;
  double *LDOS = (double *) calloc(N, sizeof(double));

  
  for (i=0, step=0; step<steps; step++) {
    // step*dE < E < (step+1)*dE
    for (j=0; j<N; j++)
      LDOS[j] = 0.0;
    while (evals[i] < (step+1)*dE) {
      for (j=0; j<N; j++)
	LDOS[j] += pow(evecs[i*N+j], 2.0);
      i++;
    } // while
    for (j=0; j<N; j++) {
      x = -L/2.0 + j*dx;
      fprintf(fp, "%g %g %g\n", x, (step+0.5)*dE*Eunit, LDOS[j]);
    } // for j
    fprintf(fp, "\n");
  } // for i


  fclose(fp);
  free(v);
  free(diag);
  free(subd);
  free(evals);
  free(evecs);
  free(ifail);
  return(0);
}
