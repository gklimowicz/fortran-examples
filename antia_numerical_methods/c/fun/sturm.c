/*	To find specified eigenvalues of a real symmetric tridiagonal 
	matrix using bisection on Sturm sequence

	E : (input) Array of length N containing the off-diagonal elements
		of the tridiagonal matrix, E[i+1]=A(i,i+1)=A(i+1,i)
	D : (input) Array of length N containing the diagonal elements
		of the tridiagonal matrix, D[i]=A(i,i)
	N : (input) Order of the matrix
	M1 : (input) Serial number of lowest eigenvalue to be determined.
		The eigenvalues are sorted in increasing order
	M2 : (input) Serial number of highest eigenvalue to be determined.
		All eigenvalues from M1 to M2 are determined
	EL : (output) Array of length M2+1 containing the calculated
		lower limit on eigenvalues
	EU : (output) Array of length M2+1 containing the calculated
		upper limit on eigenvalues
		The ith eigenvalue is located in interval (EL[i],EU[i])
	NUM : (output) Number of times Sturm sequence was evaluated to locate
		the eigenvalues.
	REPS : (input) Relative accuracy to which eigenvalues are located by bisection
		
	Error status is returned by the value of the function STURM.
		0 value implies successful execution
		110, implies that M1<0 or M2>N-1, in which case no calculations are done

	Required functions : None
*/

#include <math.h>
#include <stdlib.h>

int sturm(double e[], double d[], int n, int m1, int m2, double el[],
	double eu[], int *num, double reps)

{
	int i,k,m,nl,nu,n1,n2,ifl,ns, nbis=70;
	double emin,emax,r,e1,q1,q2, xuf=-1.0e29, eps=1.e-30;
	double *wk;

	if(m1>m2) return 0;
	if(m1<0 || m2>=n) return 110;

	e[0]=0.0;
/*	Finding bounds on eigenvalues using Gerschgorin's theorem */
	emin=d[n-1]-fabs(e[n-1]);
	emax=d[n-1]+fabs(e[n-1]);
	wk=(double *) calloc((size_t) n, sizeof(double));
	wk[0]=0.0;
	for(i=n-2; i>=0; --i) {
		r=fabs(e[i])+fabs(e[i+1]);
		if(d[i]+r>emax) emax=d[i]+r;
		if(d[i]-r<emin) emin=d[i]-r;
		wk[i+1]=e[i+1]*e[i+1];
	}

/*	Initialise the limits to undefined values */
	for(i=m1; i<=m2; ++i) {el[i]=xuf; eu[i]=xuf;}
	nl=0; nu=n;
	*num=0;

/*	Loop for each eigenvalue */
	for(m=m1; m<=m2; ++m) {
/*	If the lower bound is undefined, then use EMIN */
		if(el[m]==xuf) {el[m]=emin; n1=nl;}

		if(eu[m]==xuf) {
/*	If upper bound is undefined, use the bound for some higher eigenvalue
	and if none exists, then use EMAX */
			for(i=m+1; i<=m2; ++i) {
				if(eu[i] != xuf) {
					eu[m]=eu[i];
					n2=i;
					break;
				}
			}
			if(eu[m]==xuf) {eu[m]=emax; n2=nu;}
		}

		ifl=0;
/*	Loop for bisection */
		for(i=1; i<=nbis; ++i) {
			e1=(el[m]+eu[m])/2.;
			if(e1==el[m] || e1==eu[m]) break;
			*num=(*num)+1;

/*	Count the number of sign changes in the Sturm sequence */
			ns=0;
			q1=d[0]-e1;
			if(q1<0.0) ns=1;
			for(k=1; k<n; ++k) {
				q2=d[k]-e1;
				if(q1 != 0.0) q2=q2-wk[k]/q1;
				else q2=q2-fabs(e[k])/eps;
				if(q2<0.0) ns=ns+1;
				q1=q2;
			}

/*	If the bounds are two consecutive real number on the machine, then quit */
			if(e1 == el[m] || e1==eu[m]) break;
/*	Update the bounds */
			if(ns>=m1+1 && ns<=m2+1) {
				if(e1<eu[ns-1]) eu[ns-1]=e1;
			}
			if(ns>=m1 && ns<m2+1) {
				if(e1>el[ns]) el[ns]=e1;
			}
			if(ns<m+1) {
				if(ns>n1) n1=ns;
				if(e1>el[m]) el[m]=e1;
			}
			else {
				if(ns<n2) n2=ns;
				if(e1<eu[m]) eu[m]=e1;
			}

			if(n1==m && n2==m+1) {
/*	The eigenvalue is isolated */
				if(ifl>3 && fabs(eu[m]-el[m])) goto next;
				if(m==m1) ifl=ifl+1;
				else if(el[m] != eu[m-1]) ifl=ifl+1;
			}
		}

/*	If the eigenvalue cannot be isolated, then set the same bounds
	for all of them */
		for(k=m+1; k<n2; ++k) {el[k]=el[m]; eu[k]=eu[m];}

next:	n1=m+1; if(m+2>n2) n2=m+2;
	}
	free(wk);
	return 0;
}
