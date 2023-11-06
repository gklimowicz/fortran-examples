/*	Real eigenvalues and left eigenvector of a real matrix using inverse iteration

	A : (input) Array of length IA*M containing the matrix elements
	M : (input) Order of the matrix
	IA : (input) The second dimension of A as declared in the calling function, IA>=M
	P : (input/output) Initial value of the shift. This will be modified
		by the function if IFLG>0
	U : (input/output) Array of length M, which should specify the
		initial approximation to eigenvector. After execution it
		will contain the calculated left eigenvector.
	IFLG : (input) Integer variable to specify the type of iteration required
		If IFLG=0 the shift P is kept fixed
		If IFLG=1 the shift P is varied using Rayleigh quotient
		If IFLG=2 the shift P is varied using max(V_s+1)
	EI : (output) Estimated eigenvalue using simple inverse iteration
	ERC : (output) Estimated eigenvalue using Rayleigh quotient
	REPS : (input) Required absolute accuracy. Iteration is terminated
		when all components of eigenvector and the eigenvalue have
		converged to REPS.
	NIT : (input/output) Number of iterations required. If it is
		zero or negative NIT is set to NIT0=100
		
	Error status is returned by the value of the function INVIT_L.
		0 value implies successful execution
		106 implies that M<=1 or M>IA, in which case no
			calculations are done
		141 implies that vector is zero at some stage and
			calculations are aborted
		142 implies that inverse iteration has failed to converge

	Required functions : GAUELM
*/

#include <math.h>
#include <stdlib.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);

int invit_l(double *a, int m, int ia, double *p, double u[], int iflg,
	double *ei, double *erc, double reps, int *nit)

{
	int i,j,k,ifl,km,ier, nit0=100;
	double epi,r1,ukm,s1,s2,det;
	int *iwk;
	double *wk,*v;

	if(m<=1 || m> ia) return 106;

	iwk=(int *) calloc((size_t) m, sizeof(int));
	v=(double *) calloc((size_t) m, sizeof(double));
	wk=(double *) calloc((size_t) (m*m), sizeof(double));

/*	Copy the matrix to WK and apply the shift */
	for(i=0; i<m; ++i) {
		v[i]=u[i];
		for(j=0; j<m; ++j) wk[j+i*m]=a[i+j*m];
		wk[i+i*m]=a[i+i*m]-(*p);
	}

	ifl=1;
	if(*nit<=0) *nit=nit0;
/*	Perform Gaussian elimination on A-pI */
	ier=gauelm(m,1,wk,u,&det,iwk,m,&ifl);
	if(ier>0) {free(wk); free(v); free(iwk); return ier;}

	epi=0.0;
/*	Loop for inverse iteration */
	for(j=1; j<=(*nit); ++j) {
		ier=gauelm(m,1,wk,u,&det,iwk,m,&ifl);
		if(ier>0) {free(wk); free(v); free(iwk); return ier;}

/*	Normalising the vector U */
		r1=0.0; km=0;
		for(k=0; k<m; ++k) {
			if(r1<fabs(u[k])) {r1=fabs(u[k]); km=k;}
		}
		ukm=u[km];
/*	If the vector is zero, then quit */
		if(ukm==0.0) {free(wk); free(v); free(iwk); return 141;}

/*	The eigenvalue */
		*ei=v[km]/ukm+(*p);
		s1=0.0; s2=0.0;
/*	Calculating the Rayleigh quotient */
		for(k=0; k<m; ++k) {
			s1=s1+u[k]*v[k];
			s2=s2+fabs(u[k])*fabs(u[k]);
			u[k]=u[k]/ukm;
		}
		*erc=(*p)+s1/s2;

/*	Convergence check */
		r1=fabs(*ei-epi);
		for(i=0; i<m; ++i) {
			if(r1<fabs(v[i]-u[i])) r1=fabs(v[i]-u[i]);
			v[i]=u[i];
		}
		if(fabs(r1)<reps) {free(wk); free(v); free(iwk); return 0;}
		epi=(*ei);

		if(iflg>=1) {
/*	Update the shift */
			*p=(*erc);
			if(iflg==2) *p=(*ei);
/*	Setting up the new matrix A-pI */
			for(i=0; i<m; ++i) {
				for(k=0; k<m; ++k) wk[k+i*m]=a[i+k*m];
				wk[i+i*m]=a[i+i*m]-(*p);
			}
			ifl=1;
			ier=gauelm(m,1,wk,u,&det,iwk,m,&ifl);
			if(ier>0) {free(wk); free(v); free(iwk); return ier;}
		}
	}

	free(wk); free(v); free(iwk);
	return 142;
}
