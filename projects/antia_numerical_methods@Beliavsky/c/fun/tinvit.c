/*	To find a specified eigenvalue and eigenvector of a real symmetric
	tridiagonal matrix using inverse iteration

	E : (input) Array of length N containing the off-diagonal elements
		of the tridiagonal matrix, E[i+1]=A(i,i+1)=A(i+1,i)
	D : (input) Array of length N containing the diagonal elements
		of the tridiagonal matrix, D[i]=A(i,i)
	N : (input) Order of the matrix
	EL : (input) lower limit on eigenvalue
	EU : (input) upper limit on eigenvalue
	EI : (output) The calculated eigenvalue
	EV : (input/output) Array of length N containing the calculated
		eigenvector. If IFLG!=0 then EV should contain the previous
		eigenvector, so that initial guess to EV is orthogonal to it
	REPS : (input) Relative accuracy to which eigenvalues or eigenvectors
		are calculated
	IFLG : (input) Integer variable used as a flag to denote how initial
		approximation to eigenvector is obtained. If IFLG=0 then
		initial guess to eigenvector is generated using random numbers.
		Otherwise, initial guess is chosen to be orthogonal to EV
	NUM : (output) Number of iterations required by the function
		
	Error status is returned by the value of the function TINVIT.
		0 value implies successful execution
		144 implies that inverse iteration failed to converge

	Function RAN1(SEED) is required to generate random numbers.
	If a different function is used, the seed should be changed
	appropriately.

	Required functions : RAN1
*/

#include <math.h>
#include <stdlib.h>

double ran1(double *seed);

int tinvit(double e[], double d[], int n, double el, double eu, double *ei,
	double ev[], double reps, int iflg, int *num)

{
	int i,j,k,km, nit=100;
	double s,evs,el1,eu1,t,r,r1,s1,s2,dei,ukm,dvs,rep,hi;
	double *u, *b;
	static double seed=12345;   /*	The seed for random number generator */

	e[0]=0.0;
	u=(double *) calloc((size_t) n,sizeof(double));
	b=(double *) calloc((size_t) (4*n),sizeof(double));

	if(iflg==0) {
/*	Choose a random initial guess for eigenvector */
		for(i=0; i<n; ++i) {
			u[i]=ran1(&seed);
			ev[i]=u[i];
		}
	}

	else {
/*	Choose the initial vector to be orthogonal to EV */
		s=0.0; evs=0.0;
		for(i=0; i<n; ++i) {
			u[i]=ran1(&seed);
			evs=evs+ev[i]*ev[i];
			s=s+u[i]*ev[i];
		}
		s=s/evs;
		for(i=0; i<n; ++i) {
			ev[i]=u[i]-ev[i]*s;
			u[i]=ev[i];
		}
	}

	*ei=(el+eu)/2.;
/*	Stretch the bounds to allow for roundoff errors */
	el1=el-reps*fabs(el);
	eu1=eu+reps*fabs(eu);

/*	Loop for inverse iteration */
	for(j=1; j<=nit; ++j) {
		*num=j;

/*	Set up the matrix A-pI */
		for(i=0; i<n; ++i) {
			b[i*4]=e[i];
			b[1+i*4]=d[i]-(*ei);
			if(i<n) b[2+i*4]=e[i+1];
			b[3+i*4]=0.0;
		}
		b[4*n-2]=0.0;

/*	Gaussian elimination with partial pivoting */
		for(k=0; k<n-1; ++k) {
			if(fabs(b[1+k*4])<fabs(b[4*k+4])) {
				for(i=0; i<3; ++i) {
					t=b[i+4*k+4]; b[i+4*k+4]=b[i+1+4*k]; b[i+1+4*k]=t;
				}
				t=ev[k+1]; ev[k+1]=ev[k]; ev[k]=t;
			}

/*	If a pivot is zero make it nonzero */
			if(b[1+4*k]==0.0) b[1+4*k]=reps;
			r=-b[4*k+4]/b[1+4*k];
			b[4*k+5]=b[4*k+5]+r*b[2+4*k];
			b[4*k+6]=b[4*k+6]+r*b[3+4*k];
			ev[k+1]=ev[k+1]+r*ev[k];
		}
		if(b[4*n-3]==0.0) b[4*n-3]=reps;

/*	Back-substitution */
		ev[n-1]=ev[n-1]/b[4*n-3];
		ev[n-2]=(ev[n-2]-ev[n-1]*b[4*n-6])/b[4*n-7];
		for(k=n-3; k>=0; --k) ev[k]=(ev[k]-ev[k+1]*b[2+4*k]-ev[k+2]*b[3+4*k])/b[1+4*k];

/*	Calculating the Rayleigh quotient (DEI) */
		r1=0.0; km=0; s1=0.0; s2=0.0;
		for(k=0; k<n; ++k) {
			s1=s1+ev[k]*ev[k];
			s2=s2+ev[k]*u[k];
			if(r1<fabs(ev[k])) {r1=fabs(ev[k]); km=k;}
		}

		dei=s2/s1;
		ukm=ev[km];
		dvs=0.0;
/*	Normalising the eigenvector */
		for(k=0; k<n; ++k) {
			ev[k]=ev[k]/ukm;
			dvs=dvs+fabs(ev[k]-u[k]);
		}

/*	The convergence check */
		rep=reps*fabs(*ei);
		if(fabs(dei)<rep || r1*rep>1.0 || dvs<reps) {free(b); free(u); return 0;}
		for(i=0; i<n; ++i) u[i]=ev[i];
		if(*ei+dei>el1 && (*ei)+dei<eu1) *ei=(*ei)+dei;
		else {
/*	If the new shift is outside the specified bounds, then modify it */
			hi=fabs(*ei-el); if(hi>fabs(*ei-eu)) hi=fabs(*ei-eu);
			s1=hi*0.5; if(dei<0.0) s1=-s1;
			*ei=(*ei)+s1;
		}
	}

	free(b); free(u);
	return 144;
}
