/*	Inverse Laplace transform */

#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979324

int lapinv(int n, double t[], double f[], void cfs(double * , double * ),
	double alpha, double *reps);
void cf(double s[2], double f[2]);
void cdiv(double c1[2], double c2[2], double cr[2]);

main()
{
	int i,i1,j,n, id, iflg, ier,np;
	double hh,x,a,f[2000],t[2000],reps,tmax;

/*	Example 10.7 */

	np=100; iflg=1; reps=1.e-6;

/*	The program chooses n uniformly spaced points in interval [0, tmax]
	to evaluate the result */

	for(i1=0; i1<99; ++i1) {
		printf("type n=no. of points,  a=exponential order,  tmax \n");
		printf("                   (quits when n<=0)\n");
		scanf(" %d %le %le", &n, &a, &tmax);
		if(n<=0) return 0;

		hh=tmax/(n-1);
		for(i=0; i<n; ++i) {
			t[i]=i*hh;
		}

		i=lapinv(n,t,f,cf,a,&reps);
		printf(" ier = %d  n = %d  exponential order = %e\n",i,n,a);
		for(i=0; i<n; ++i) printf(" f(%e) = %e \n",t[i],f[i]);

	}
	return;
}

void cf(double cs[2], double f[2])

{
	double c1[2],c2[2],c3[2],c4[2];

	c1[0]=2.0; c1[1]=0.0; cdiv(c1,cs,c2);
	c1[0]=1.0; c4[0]=cs[0]+1.0; c4[1]=cs[1]; cdiv(c1,c4,c3);
	f[0]=c2[0]-c3[0]; f[1]=c2[1]-c3[1];
}




/*	To calculate inverse Laplace transform

	N : (input) Number of points at which inverse transform needs to
		be evaluated
	T : (input) Array of length N, containing the points at
		which inverse transform needs to be calculated. These
		elements can be in any order but the last element must
		be the maximum.
	F : (output) Array of length N, containing the calculated
		function values. F[I] is the inverse Laplace transform at T[I]
	CFS : (input) Name of the function to calculate the
		Laplace transform, which is to be inverted
	ALPHA : (input) The estimated exponential order of the function
	REPS : (input/output) The required relative accuracy
		
	Error status is returned by the value of the function LAPINV.
		0 value implies successful execution
		61 implies that epsilon-algorithm failed for at least
			one value of T[I]
		62 implies that denominator is zero for some T[I],
			in this case the last value is accepted

	Function CFS(S,FS) must be supplied by the user. Here both S and
		and FS are arrays of length 2 containing the value of complex
		variables. FS is f(S).

	Required functions : CFS
*/

#include <math.h>

#define PI 3.14159265358979324

int lapinv(int n, double t[], double f[], void cfs(double * , double * ),
	double alpha, double *reps)

{
	int i,j,nf,it,nd,nt,i1,i2,ii,ier, nmin=3, nmax=100;
	double tn,a,pit,pi0,f1,s1,t1,t2,dif,dif1,ri,den;
	double cf[2],cs[2],cwk[101][2],et[2][100];

	tn=t[n-1];
	if(*reps<=0.0) *reps=1.e-6;
	a=alpha-0.5*log(0.5*(*reps))/tn;
	pit=PI/tn;
	cs[0]=a; cs[1]=0.0;
	cfs(cs,&cwk[0][0]);
	nf=0;
	ier=0;

/*	Loop over the N points */
	for(it=0; it<n; ++it) {
		pi0=pit*t[it];
		f1=exp(a*t[it])/tn;
		s1=0.5*cwk[0][0];
		nd=1; nt=2;
		et[0][0]=0.0;
		i1=0; i2=1;

		for(i=1; i<nmax; ++i) {
			for(j=nd; j<=nt; ++j) {
				cs[0]=a; cs[1]=j*pit;
				if(j>nf) {
/*	Evaluate the function */
					cfs(cs,cf); nf=nf+1;
					cwk[nf][0]=cf[0]; cwk[nf][1]=cf[1];
				}
				else {
/*	else use the old value of function */
					cf[0]=cwk[j][0]; cf[1]=cwk[j][1];
				}

				t2=cf[0];
				t1=t2*cos(j*pi0)-cf[1]*sin(j*pi0);
				s1=s1+t1;
			}
/*	The second column of epsilon-table */
			et[i2][1]=s1; et[i2][0]=0.0;
			if(i==2) {dif=fabs(et[i2][1]-et[i1][1]); ri=et[i2][1];}

/*	The epsilon-algorithm (ET[J][I]=epsilon(j-i+2,i-2)) */
			for(j=2; j<i+1; ++j) {
				den=et[i2][j-1]-et[i1][j-1];
				if(den != 0.0) et[i2][j]=et[i1][j-2]+1.0/den;
				else {ier=62; et[i2][j]=et[i1][j-2];}
			}

			if(i>2) {
/*	Perform convergence check */
				for(j=1; j<i-1; j += 2) {
					dif1=fabs(et[i2][j]-et[i1][j]);
					if(dif1<dif) {dif=dif1; ri=et[i2][j];}
				}
			}

			nd=nt+1;
			if(100.0*dif<(*reps)*fabs(ri) && i>nmin) break;

			ii=i1; i1=i2; i2=ii;
			nt=nd;
		}

		if(i>=nmax) ier=61;
		f[it]=f1*ri;
	}
	return ier;
}




/* Utility functions for complex arithmetic, complex numbers are
	represented as double array of length 2, with first element as
	real part and second element as imaginary part.

*/

/*	Complex division

	C1 : (input) Array of length 2 containing the numerator 
	C2 : (input) Array of length 2 containing the denominator
	CR : (output) Array of length 2 containing the result (C1/C2)
		CR can be same as  C1 or C2

	Required functions : None

*/

#include <math.h>

void cdiv(double c1[2], double c2[2], double cr[2])

{
	double r,r1,den;

	if(fabs(c2[0])>fabs(c2[1])) {
		r=c2[1]/c2[0];
		den=c2[0]+c2[1]*r;
/*	To avoid overwriting on c1 if c1 and cr are same */
		r1=(c1[0]+c1[1]*r)/den;
		cr[1]=(c1[1]-c1[0]*r)/den;
		cr[0]=r1;
	}
	else {
		r=c2[0]/c2[1];
		den=c2[0]*r+c2[1];
		r1=(c1[0]*r+c1[1])/den;
		cr[1]=(c1[1]*r-c1[0])/den;
		cr[0]=r1;
	}
	return;
}
