/*	Multiple integration in n dimensions */

#include <stdio.h>
#include <math.h>

double fun(int n, double x[]);
int ngauss(double a[], double b[], int n, int m[], int ind[], double (*f) (int , double * ),
	double *ri, int *num, int *maxpt);
int mulint(double a[], double b[], int n, int m[], int ind[], double (*f) (int , double * ),
	double *rint, double reps, double aeps, double *dif, int *num, int *maxpt);
int stroud(double a[], double b[], int n, int m, int ind[], double (*f) (int , double * ),
	double *ri, int *num, int *maxpt);
int strint(double a[], double b[], int n, int *m, int ind[], double (*f) (int , double * ),
	double *rint, double reps, double aeps, double *dif, int *num, int *maxpt);
double ranf(int *iseed);
int mcarlo(double a[], double b[], int n, int npt, double (*f) (int , double * ),
	double *ri, double reps, double aeps, double *err, int *np);
int equids(double a[], double b[], int n, int npt, double (*f) (int , double * ),
	double *s1, double *s2, double reps, double aeps, double *dif, int *np);


main()
{
	int i,i1,j,n,nuse, m[20],ind[20], iflg,m0, ier,np,nmax;
	double a[20], b[20], rint, reps, aeps, df, s2;

/*	Example 6.15   (I3)  */

	aeps=1.e-8; reps=1.e-7;
	nmax=1000000;
	for(i=0; i<20; ++i) {a[i]=0.0; b[i]=1.0;}
	for(i1=0; i1<99; ++i1) {
		printf("type n = No. of dimensions,  m0 = order of formula (stroud) \n");
		printf("                         (quits when n<=0)\n");
		scanf(" %d %d", &n, &m0);
		if(n<=0) return 0;

		for(i=0; i<n; ++i) {m[i]=0; ind[i]=0;}
		i=mulint(a,b,n,m,ind,fun,&rint,reps,aeps,&df,&nuse,&nmax);
		printf("mulint:  ier = %d  n = %d  no. of function evaluations = %d \n", i,n,nuse);
		printf(" integral = %e   estimated error = %e \n", rint,df);

		for(i=0; i<n; ++i) {ind[i]=0;}
		i=strint(a,b,n,&m0,ind,fun,&rint,reps,aeps,&df,&nuse,&nmax);
		printf("strint:  ier = %d  no. of function evaluations = %d  order = %d\n", i,nuse,m0);
		printf(" integral = %e   estimated error = %e \n", rint,df);

		i=mcarlo(a,b,n,nmax,fun,&rint,reps,aeps,&df,&nuse);
		printf("mcarlo:   ier = %d  no. of function evaluations = %d \n", i,nuse);
		printf(" integral = %e   estimated error = %e \n", rint,df);

		i=equids(a,b,n,nmax,fun,&rint,&s2,reps,aeps,&df,&nuse);
		printf("equids:   ier = %d  no. of function evaluations = %d\n", i,nuse);
		printf(" integral (s1,s2) = %e  %e   estimated error = %e \n", rint,s2,df);

	}
	return;
}

double fun(int n, double x[])

{
	int i;
	double xa;

	xa=0.0;
	for(i=0; i<n; ++i) xa=xa+x[i]*x[i];
	if(n-xa > 0.0) return 1./sqrt(n-xa);
	else return 0.0;
}




/*	Multiple integration over a hyper-rectangle in n-dimensions
	using equidistributed sequences

	A : (input) Array of length N containing the lower limit
		along each dimension
	B : (input) Array of length N containing the upper limit
		along each dimension
	N : (input) The number of dimensions
	NPT : (input) Maximum number of function evaluations to be used
	F : (input) Name of the function to calculate the integrand
		Function F(N,X) should calculate the integrand, where N is the
		number of dimensions and X is an array of length N containing
		the coordinates of the point where integrand is to be calculated
	S1 : (output) The calculated value of the integral
	S2 : (output) Another approximation to the value of the integral
		For smooth functions S2 is expected to be better approximation
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(S2))
	DIF : (output) estimated (absolute) error achieved by the function
	NP : (output) Number of function evaluations used
		
	Error status is returned by the value of the function EQUIDS.
		0 value implies successful execution
		39 implies specified accuracy was not achieved in
			which case DIF will contain the estimated accuracy
		312 implies N<1 or N>21 and no calculations are done

	Function F(N,X) must be supplied by the user
	
	Required functions :  F
*/

#include <math.h>


int equids(double a[], double b[], int n, int npt, double (*f) (int , double * ),
	double *s1, double *s2, double reps, double aeps, double *dif, int *np)

{
	int i,j,npt1,ia, nchk=100, nmax=21;
	double ri,hh,ri1,dif1,diff,ss1,a1,r1,h[21],xa[21],wt[21];
/*	The first 21 prime numbers */
	double at[21] ={2.,3.,5.,7.,11.,13.,17.,19.,23.,29.,31.,37.,41.,
             43.,47.,53.,59.,61.,67.,71.,73.};

	*s1=f(n,a);
	*s2=(*s1);
	*np=1;
	if(n>nmax || n<1) return 312;

	hh=1.0;
	for(i=0; i<n; ++i) {
		h[i]=b[i]-a[i];
		hh=hh*h[i];
/*	The irrational numbers for generating equidistributed sequences */
		wt[i]=sqrt(at[i]);
	}

	ri=0.0; ri1=0.0; *dif=0.0;
	npt1=nchk;
	for(i=1; i<=npt; ++i) {
/*	Generate the abscissas using equidistributed sequences */
		for(j=0; j<n; ++j) {
			a1=i*wt[j];
			ia=a1+0.5; a1=2.*fabs(a1-ia)*h[j];
			xa[j]=a[j]+a1;
		}
/*	Accumulate the sum */
		*s1=(*s1)+2.*(f(n,xa)-ri);
		*s2=(*s2)+(*s1);

		if(i-nchk*(i/nchk) == 0) {
/*	To control the roundoff error form partial sums */
			ss1=ri+(*s1)/(2*i+1);
			diff=(*s2)/((i+1.)*(i+1.));
			*s2=0.0;
			*s1=(*s1)-(2*i+1)*diff;
/*	The new approximation to the average value of function */
			ri=ri+diff;

			if(i==npt1) {
/*	Check for convergence */
				dif1=*dif;
				*dif=fabs(ri-ri1);
				ri1=ri;
				r1=reps*fabs(ri); if(aeps/hh > r1) r1=aeps/hh;
				if((*dif+dif1 < r1) && i > 5*nchk) {
					*s1=ss1*hh;
					*s2=ri*hh;
					*dif=(*dif+dif1)*hh;
					*np=i+1;
					return 0;
				}

				npt1=npt1*2;
			}
		}
	}

/*	Integral fails to converge */
	*s1=(ri+(*s1)/(2*npt+1))*hh;
	*s2=(ri+(*s2)/((npt+1.0)*(npt+1.0)))*hh;
	*dif=(*dif+dif1)*hh;
	*np=npt+1;
	return 39;
}



/*	Multiple integration over a hyper-rectangle in n-dimensions
	using Monte-Carlo technique

	A : (input) Array of length N containing the lower limit along each dimension
	B : (input) Array of length N containing the upper limit along each dimension
	N : (input) The number of dimensions
	NPT : (input) Maximum number of function evaluations to be used
	F : (input) Name of the function to calculate the integrand
		Function(N,X) should calculate the integrand, where N is the
		number of dimensions and X is an array of length N containing
		the coordinates of the point where integrand is to be calculated
	RI : (output) The calculated value of the integral
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RI))
	ERR : (output) estimated (absolute) error achieved by the function
		It is 2.576 times the estimated standard deviation
		divided by SQRT(NP).
	NP : (output) Number of function evaluations used
		
	Error status is returned by the value of the function MCARLO.
		0 value implies successful execution
		39 implies specified accuracy was not achieved in
			which case ERR will contain the estimated accuracy
		311 implies N<1, and no calculations are done

	Function F(N,X) should be supplied by the user
	
	Required functions : RANF, F
*/

#include <math.h>
#include <stdlib.h>

double ranf(int *iseed);
double ran(double *seed);


int mcarlo(double a[], double b[], int n, int npt, double (*f) (int , double * ),
	double *ri, double reps, double aeps, double *err, int *np)

{
	int i,j,iseed,npt1, nchk=100;
	double ri1,var,var1,hh,f1,r1;
	double *h, *xa;

	if(n<1) return 311;
	*ri=0.0; *np=0;
	h=(double *) calloc((size_t) n, sizeof(double));
	xa=(double *) calloc((size_t) n, sizeof(double));

	hh=1.0;
	for(i=0; i<n; ++i) {h[i]=b[i]-a[i]; hh=hh*h[i];}

	ri1=0.0; var1=0.0;

/*	Seed for random number generator, should be changed if another routine is used */
	iseed=-12345;
	npt1=nchk;

	for(i=1; i<=npt; ++i) {
/*	Generating the abscissas */
		for(j=0; j<n; ++j) xa[j]=a[j]+h[j]*ranf(&iseed);
		f1=f(n,xa);
		ri1=ri1+f1;
		var1=var1+f1*f1;

		if(i == npt1) {
/*	Compute intermediate sums to check for convergence */
			*ri=ri1/i;
			var=var1/i-(*ri)*(*ri);
			if(var<0.0) var=0.0;
			*err=2.576*hh*sqrt(var/npt1);
			*ri=(*ri)*hh;
			*np=i;
			r1=reps*fabs(*ri); if(aeps>r1) r1=aeps;
			if(*err < r1) {free(xa); free(h); return 0;}
			npt1=2*npt1;
		}
	}

/*	Integral fails to converge */
	*ri=ri1/npt;
	var=var1/npt-(*ri)*(*ri);
	*err=2.576*hh*sqrt(var/npt);
	*ri=(*ri)*hh;
	*np=npt;
	r1=reps*fabs(*ri); if(aeps>r1) r1=aeps;
	if(*err > r1) {free(xa); free(h); return 39;}
	free(xa); free(h);
	return 0;
}




/*	Multiple integration over a hyper-rectangle in n-dimensions
	using product Gauss-Legendre formulas

	A : (input) Array of length N containing the lower limit along each dimension
	B : (input) Array of length N containing the upper limit along each dimension
	N : (input) The number of dimensions
	M : (input/output) Integer array of length N specifying the formula
		to be used along each dimension. M[J]-point formula will
		be used along Jth dimension, M[J] should be 2,4,8,16 or 32
		otherwise it will be set to a default value of 2
	IND : (input/output) Integer array of length N specifying the number
		of subintervals to be used along each dimension. IND[J]>0
		otherwise it will be set to a default value of 1
	F : (input) Name of the function to calculate the integrand
		Function(N,X) should calculate the integrand, where N is the
		number of dimensions and X is an array of length N containing
		the coordinates of the point where integrand is to be calculated
	RINT : (output) The calculated value of the integral
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RINT))
	DIF : (output) estimated (absolute) error achieved by the function
	NUM : (output) Number of function evaluations used
	MAXPT : (input/output) Maximum number of function evaluations permitted
		If MAXPT <1 it is set to a default value of MAXPTS (=1100000)
		
	Error status is returned by the value of the function MULINT.
		0 value implies successful execution
		39 implies specified accuracy was not achieved in
			which case DIF will contain the estimated accuracy
		305 implies N<1 and no calculations are done
		307 implies that number of points exceeded MAXPT in
			first attempt and no approximation of RINT is calculated

	Function F(N,X) must be supplied by the user
	
	Required functions : NGAUSS, F
*/

#include <math.h>

int ngauss(double a[], double b[], int n, int m[], int ind[], double (*f) (int , double * ),
	double *ri, int *num, int *maxpt);

int mulint(double a[], double b[], int n, int m[], int ind[], double (*f) (int , double * ),
	double *rint, double reps, double aeps, double *dif, int *num, int *maxpt)

{
	int i,j,i1,m1,no,qc,ier, maxpts=1100000;
	double rint1,r1;

/*	Set M[I] and IND[I] to default values if they are unacceptable */
	for(i=0; i<n; ++i) {
		if((m[i]!=2) && (m[i]!=4) && (m[i]!=8) && (m[i]!=16) && (m[i]!=32)) m[i]=2;
		if(ind[i] <1) ind[i]=1;
	}
	if(*maxpt <1) *maxpt=maxpts;
	*num=0;

/*	Evaluate the integral */
	ier=ngauss(a,b,n,m,ind,f,rint,&no,maxpt);
	*num=(*num)+no;
	if(ier > 100) return ier;

/*	Iteration to check and improve the accuracy of integral */
	for(i=0; i<10; ++i) {
		qc=1;
		*dif=0.0;

/*	Check for convergence along each dimension by doubling the
	number of points at which function is evaluated */
		for(j=0; j<n; ++j) {
			if(*num > *maxpt) return 39;
			m1=m[j]; i1=ind[j];
/*	If M[J]<32 then double the order of formula */
			if(m[j]<32) m[j]=2*m[j];

/*	otherwise double the number of subintervals */
			else ind[j]=2*ind[j];

			ier=ngauss(a,b,n,m,ind,f,&rint1,&no,maxpt);
			*num=(*num)+no;
			*dif=(*dif)+fabs(*rint-rint1);
			if(ier > 100) return 39;

			r1=reps*fabs(rint1); if(aeps>r1) r1=aeps;
			if(fabs(rint1-*rint) < r1) {m[j]=m1; ind[j]=i1;}
			else {*rint=rint1; qc=0;}
		}

/*	If satisfactory accuracy is achieved for all dimensions then return */
		if(qc==1) return 0;

/*	If the number of function evaluations exceeds maxpt then quit */
		if(*num> (*maxpt)) return 39;
	}
	return 39;
}



/*	Multiple integration over a hyper-rectangle in n-dimensions
	using product Gauss-Legendre formulas with given number of points

	A : (input) Array of length N containing the lower limit
		along each dimension
	B : (input) Array of length N containing the upper limit
		along each dimension
	N : (input) The number of dimensions, N>0
	M : (input) Integer array of length N specifying the formula
		to be used along each dimension. M[J]-point formula will
		be used along Jth dimension, M[J] should be 2,4,8,16 or 32
		otherwise no calculations are done
	IND : (input) Integer array of length N specifying the number
		of subintervals to be used along each dimension. IND[J]>0
		otherwise no calculations are done
	F : (input) Name of the function to calculate the integrand
		Function(N,X) should calculate the integrand, where N is the
		number of dimensions and X is an array of length N containing
		the coordinates of the point where integrand is to be calculated
	RI : (output) The calculated value of the integral
	NUM : (output) Number of function evaluations used
	MAXPT : (input/output) Maximum number of function evaluations permitted
		If MAXPT <1 it is set to a default value of MAXPTS (=1100000)
		
	Error status is returned by the value of the function NGAUSS.
		0 value implies successful execution
		305 implies N<1, in which case no calculations are done
		306 implies M[J] is not 2,4,8,16 or 32 or IND[J]<1 for some J
			in which case no calculations are done
		307 implies that number of points exceeded MAXPT and
			no calculations are done
	
	Required functions : F
*/

#include <math.h>
#include <stdlib.h>

#define MAXPTS 1100000

int ngauss(double a[], double b[], int n, int m[], int ind[], double (*f) (int , double * ),
	double *ri, int *num, int *maxpt)

{
	int i,j,no,j1,j2,m1,m2,k,ier, maxpts=1100000;
	double x1,h1;
	int *ip, *npt;
	double *h, *xa, *wt;

/*	Weights and abscissas for Gauss-Legendre quadrature.
	For N-point formula W[K]=W[N-K-1] and X[K]=-X[N-K-1]
		For K=0,1,...,N/2-1. Hence only half points are tabulated.
	For 2-point W[0]; 4-point W[1], W[2]; 8-point W[3],...,W[6];
	16-point W[7],...,W[14]; 32-point W[15],...,W[30] are the
	weights corresponding to abscissas X[I].
*/

	double w[31] = {1.0e0,
      0.34785484513745385737e0, 0.65214515486254614263e0,
      0.10122853629037625915e0, 0.22238103445337447054e0,
      0.31370664587788728734e0, 0.36268378337836198297e0,
      0.02715245941175409485e0, 0.06225352393864789286e0,
      0.09515851168249278481e0, 0.12462897125553387205e0,
      0.14959598881657673208e0, 0.16915651939500253819e0,
      0.18260341504492358887e0, 0.18945061045506849629e0,
      0.00701861000947009660e0, 0.01627439473090567061e0,
      0.02539206530926205945e0, 0.03427386291302143310e0,
      0.04283589802222668066e0, 0.05099805926237617620e0,
      0.05868409347853554714e0, 0.06582222277636184684e0,
      0.07234579410884850623e0, 0.07819389578707030647e0,
      0.08331192422694675522e0, 0.08765209300440381114e0,
      0.09117387869576388471e0, 0.09384439908080456564e0,
      0.09563872007927485942e0, 0.09654008851472780057e0};

	double x[31] = {0.57735026918962576451e0,
      0.86113631159405257522e0, 0.33998104358485626480e0,
      0.96028985649753623168e0, 0.79666647741362673959e0,
      0.52553240991632898582e0, 0.18343464249564980494e0,
      0.98940093499164993260e0, 0.94457502307323257608e0,
      0.86563120238783174388e0, 0.75540440835500303390e0,
      0.61787624440264374845e0, 0.45801677765722738634e0,
      0.28160355077925891323e0, 0.09501250983763744019e0,
      0.99726386184948156354e0, 0.98561151154526833540e0,
      0.96476225558750643077e0, 0.93490607593773968917e0,
      0.89632115576605212397e0, 0.84936761373256997013e0,
      0.79448379596794240696e0, 0.73218211874028968039e0,
      0.66304426693021520098e0, 0.58771575724076232904e0,
      0.50689990893222939002e0, 0.42135127613063534536e0,
      0.33186860228212764978e0, 0.23928736225213707454e0,
      0.14447196158279649349e0, 0.04830766568773831623e0};

	*ri=0.0; *num=0;
	if(n<1) return 305;
	if(*maxpt < 1) *maxpt=maxpts;

/*	calculate the number of function evaluations required */
	*num=m[0]*ind[0];
	for(i=1; i<n; ++i) *num=(*num)*m[i]*ind[i];
	if(*num > (*maxpt)) return 307;

/*	Initialisation */
	ip=(int *) calloc((size_t) (n+2), sizeof(int));
	npt=(int *) calloc((size_t) (n+2), sizeof(int));
	h=(double *) calloc((size_t) (n+2), sizeof(double));
	xa=(double *) calloc((size_t) (n+2), sizeof(double));
	wt=(double *) calloc((size_t) (n+2), sizeof(double));
	ier=0;
	for(i=0; i<n; ++i) {
		ip[i]=0;
		npt[i]=m[i]*ind[i]-1;
		if((m[i]!=2) && (m[i]!=4) && (m[i]!=8) && (m[i]!=16) && (m[i]!=32)) ier=306;
		if(ind[i] < 1) ier=306;
		h[i]=(b[i]-a[i])/(2*ind[i]);
	}
	if(ier!=0) {free(wt); free(xa); free(h); free(npt); free(ip); return ier;}
	for(i=n; i<n+1; ++i) {h[i]=1.0; wt[i]=1.0;}

/*	Loop for sum over N dimensions */
	k=n-1;
	do {
		for(i=k; i>=0; --i) {
			m2=m[i]/2;
/*	The abscissas are X(NO),...,X(NO+M2-1) */
			no=m2-1;
			h1=h[i];
			j1=ip[i]/m[i];
			j2=ip[i]-j1*m[i];
/*	Use the (J2)th point in (J1)th subinterval */
			x1=a[i]+(2*j1+1)*h1;
			if(j2<m2) {
/*	For the first M2 abscissas */
				xa[i]=x1+h1*x[no+j2];
				wt[i]=w[no+j2]*wt[i+1];
			}
			else if(j2-m2 < m2) {
/*	For the next M2 abscissas */
				xa[i]=x1-h1*x[no+j2-m2];
				wt[i]=w[no+j2-m2]*wt[i+1];
			}
			else {
/*	For Gaussian formula with odd number of points use the abscissa at x=0 */
				xa[i]=x1;
				wt[i]=w[no+m2]*wt[i+1];
			}
		}

/*	Add the new point to running sum */
		*ri=(*ri)+wt[0]*f(n,xa);
		k=0;
		while(k <= n-1) {
			if(ip[k] >= npt[k]) {
				ip[k]=0;
				k=k+1;
			}
			else {
				ip[k]=ip[k]+1;
				break;
			}
		}
	} while(k<n);

/*	If all points are exhausted compute the value of integral */
	for(i=0; i<n; ++i) *ri=(*ri)*h[i];
	free(wt); free(xa); free(h); free(npt); free(ip);
	return 0;
}



/*	Multiple integration over a hyper-rectangle in n-dimensions
	using compound monomial rules

	A : (input) Array of length N containing the lower limit
		along each dimension
	B : (input) Array of length N containing the upper limit
		along each dimension
	N : (input) The number of dimensions
	M : (input/output) Integer specifying the formula to be used
		M can be 1, 3 or 5, otherwise it will be set to a default value of 3
		M=1 selects 1-point formula of degree 1
		M=3 selects 2N-point formula of degree 3 due to Stroud
		M=5 selects (2N*N+1)-point formula of degree 5
	IND : (input/output) Integer array of length N specifying the number
		of subintervals to be used along each dimension. IND[J]>0
		otherwise it will be set to a default value of 1
	F : (input) Name of the function to calculate the integrand
		Function(N,X) should calculate the integrand, where N is the
		number of dimensions and X is an array of length N containing
		the coordinates of the point where integrand is to be calculated
	RINT : (output) The calculated value of the integral
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RINT))
	DIF : (output) estimated (absolute) error achieved by the function
	NUM : (output) Number of function evaluations used by STRINT
	MAXPT : (input/output) Maximum number of function evaluations permitted
		If MAXPT <1 it is set to a default value of MAXPTS (=1000000)
		
	Error status is returned by the value of the function STRINT.
		0 value implies successful execution
		39 implies specified accuracy was not achieved in
			which case DIF will contain the estimated accuracy
		308 implies that number of points exceeded MAXPT in
			first attempt and no approximation of RINT is calculated
		309 implies N<1 or N>50 and no calculations are done

	Function F(N,X) must be supplied by the user
	
	Required functions : STROUD, F
*/

#include <math.h>

int stroud(double a[], double b[], int n, int m, int ind[], double (*f) (int , double * ),
	double *ri, int *num, int *maxpt);


int strint(double a[], double b[], int n, int *m, int ind[], double (*f) (int , double * ),
	double *rint, double reps, double aeps, double *dif, int *num, int *maxpt)

{
	int i,j,no,qc,ier,i1, maxpts=1100000;
	double rint1,r1;

/*	set M and IND[I] to default values if they are unacceptable */
	if((*m!=1) && (*m!=3) && (*m!=5)) *m=3;
	for(i=0; i<n; ++i) {
		if(ind[i]<1) ind[i]=1;
	}
	if(*maxpt <= 0) *maxpt=maxpts;
	*num=0;

/*	Evaluate the first approximation to the integral */
	ier=stroud(a,b,n,*m,ind,f,rint,&no,maxpt);
	*num= (*num)+no;
	if(ier >100) return ier;

/*	Iteration to check and improve the accuracy of integral */
	for(i=0; i<10; ++i) {
		qc=1; *dif=0.0;

		for(j=0; j<n; ++j) {
			if(*num > (*maxpt)) return 39;
			i1=ind[j];
/*	Double the number of subintervals along the Jth axis */
			ind[j]=2*ind[j];

			ier=stroud(a,b,n,*m,ind,f,&rint1,&no,maxpt);
			*num=(*num)+no;
			*dif=(*dif)+fabs(*rint-rint1);
			if(ier>100) return 39;

			r1=reps*fabs(rint1); if(aeps>r1) r1=aeps;
/*	If satisfactory accuracy is achieved then restore the old value of IND[J] */
			if(fabs(rint1-(*rint)) < r1) ind[j]=i1;
/*	else retain the new value */
			else {*rint=rint1; qc=0;}
		}

		if(qc==1) return 0;
		if(*num > (*maxpt)) return 39;
	}

	return 39;
}



/*	Multiple integration over a hyper-rectangle in n-dimensions
	using compound monomial rules with given number of points

	A : (input) Array of length N containing the lower limit
		along each dimension
	B : (input) Array of length N containing the upper limit
		along each dimension
	N : (input) The number of dimensions, N>0 and N<=NMAX (=50)
	M : (input) Integer specifying the formula to be used
		M can be 1,3 or 5, otherwise no calculations are done
		M=1 selects 1-point formula of degree 1
		M=3 selects 2N-point formula of degree 3 due to Stroud
		M=5 selects (2N*N+1)-point formula of degree 5
	IND : (input) Integer array of length N specifying the number
		of subintervals to be used along each dimension. IND[J]>0
		otherwise no calculations are done
	F : (input) Name of the function to calculate the integrand
		Function(N,X) should calculate the integrand, where N is the
		number of dimensions and X is an array of length N containing
		the coordinates of the point where integrand is to be calculated
	RI : (output) The calculated value of the integral
	NUM : (output) Number of function evaluations used by STROUD
	MAXPT : (input/output) Maximum number of function evaluations permitted
		If MAXPT <1 it is set to a default value of MAXPTS (=1000000)
		
	Error status is returned by the value of the function STROUD.
		0 value implies successful execution
		308 implies that number of points exceeded MAXPT and
			no calculations are done
		309 implies N<1 or N>=NMAX, in which case no calculations are done
		310 implies M is not 1,3 or 5 or IND[J]<1 for some J
			in which case no calculations are done

	Function F(N,X) must be provided by the user
	
	Required functions :  F
*/

#include <math.h>

#define PI 3.14159265358979324e0

int stroud(double a[], double b[], int n, int m, int ind[], double (*f) (int , double * ),
	double *ri, int *num, int *maxpt)

{
	int i,j,k,mpt,in,ier,ip[50], nmax=50, maxpts=1100000;
	double r,s1,s2,a0,a1,a2,an,xi,xi2,xi3;
	double x[50],x3[50][50],h[50],xa[50],wt[50];

	*ri=0.0; *num=0;
	if(n>nmax || n<1) return 309;

/*	Calculate the number of function evaluations required */
	if(m == 1) mpt=1;
	else if(m == 3) mpt=2*n;
	else if(m == 5) mpt=2*n*n+1;
	else return 310;

	*num=mpt;
	for(i=0; i<n; ++i) {
		if(ind[i]<1) return 310;
		*num=(*num)*ind[i];
	}
	if(*maxpt <1) *maxpt=maxpts;
	if(*num > (*maxpt)) return 308;

/*	Constants for the (2N*N+1)-point formula of degree 5 */
	xi=sqrt(0.6e0);
	a0=(25*n*n-115*n+162.0)/162.0;
	a1=(70.0-25*n)/162.0;
	a2=25.0/324.0;

/*	Abscissas for the 2N-point formula of degree 3 */
	if(m == 3) {
		xi3=sqrt(2.0/3.0);
		xi2=1.0/sqrt(3.0);
		for(i=0; i<n; ++i) {
			for(j=0; j<n-1; j=j+2) {
				an=(j+1)*(i+1)*PI/n;
				x3[j][i]=xi3*cos(an);
				x3[j+1][i]=xi3*sin(an);
			}
/*	When N is odd */
			if(n != 2*(n/2)) {
				if(i == 2*(i/2)) x3[n-1][i]=-xi2;
				else x3[n-1][i]=xi2;
			}
		}
	}

	for(i=0; i<n; ++i) {
		ip[i]=1;
		h[i]=(b[i]-a[i])/(2*ind[i]);
/*	For abscissas of (2N*N+1)-point formula of degree 5 */
		wt[i]=h[i]*xi;
	}

/*	loop for the sum over all subintervals */
	k=n-1;
	do {
		for(in=k; in>=0; --in) xa[in]=a[in]+(2*ip[in]-1)*h[in];

/*	Generalised midpoint rule */
		if(m == 1) r=f(n,xa);

		else if(m==3) {
/*	Stroud's 2N-point rule of degree 3 */
			r=0.0;
			for(i=0; i<n; ++i) {
				for(j=0; j<n; ++j) x[j]=xa[j]+x3[j][i]*h[j];
				r=r+f(n,x);
				for(j=0; j<n; ++j) x[j]=xa[j]-x3[j][i]*h[j];
				r=r+f(n,x);
			}
			r=r/(2*n);
		}

		else if(m == 5) {
/*	(2N*N+1)-point rule of degree 5 */
			r=a0*f(n,xa);
			s1=0.0; s2=0.0;
			for(i=0; i<n; ++i) x[i]=xa[i];
			for(i=0; i<n; ++i) {
				x[i]=x[i]+wt[i];
				s1=s1+f(n,x);
				x[i]=xa[i]-wt[i];
				s1=s1+f(n,x);
				x[i]=xa[i];
				for(j=i+1; j<n; ++j) {
					x[i]=xa[i]+wt[i];
					x[j]=xa[j]+wt[j];
					s2=s2+f(n,x);
					x[j]=xa[j]-wt[j];
					s2=s2+f(n,x);
					x[i]=xa[i]-wt[i];
					s2=s2+f(n,x);
					x[j]=xa[j]+wt[j];
					s2=s2+f(n,x);
					x[j]=xa[j];
					x[i]=xa[i];
				}
			}
			r=r+a1*s1+a2*s2;
		}

		*ri=(*ri)+r;
		k=0;
		while(k <= n-1) {
			if(ip[k]>= ind[k]) {
				ip[k]=1;
				k=k+1;
			}
			else {
				ip[k]=ip[k]+1;
				break;
			}
		}
	} while (k<n);

/*	If all directions are exhausted, compute the value of integral */
	for(i=0; i<n; ++i) *ri=2.*(*ri)*h[i];
	return 0;
}




/*	To generate uniformly distributed random numbers in interval (0,1)

	ISEED : (input/output) is an integer value used as the seed
		It should be initialised to negative value before first call
		and should not be modified between successive calls.

	Required functions : None
*/

#include <math.h>


double ranf(int *iseed)

{
	int m1=714025, ia1=1366, ic1=150889;
	int m2=214326, ia2=3613, ic2=45289;
	int m3=139968, ia3=3877, ic3=29573;
	int i,j, ish=43;
	double r1;

	static double rm1,rm2,ran[43];
	static int iflg=0, is1,is2,is3;
	
/*	Initialise on first call or when ISEED<0 */
	if(*iseed < 0 || iflg == 0) {
		iflg=1;
		rm1=1.0/m1;
		rm2=1.0/m2;

/*	Seeds for the three random number generators */
		is1=-(*iseed); is1=is1-m1*(is1/m1);
		is2=ia1*is1+ic1; is2=is2-m1*(is2/m1);
		is3=ia2*is2+ic2; is3=is3-m2*(is3/m2);
		*iseed=1;

/*	Store ISH random numbers in the array RAN */
		for(j=0; j<ish; ++j) {
			is1=ia1*is1+ic1; is1=is1-m1*(is1/m1);
			is2=ia2*is2+ic2; is2=is2-m2*(is2/m2);
			ran[j]=(is1+is2*rm2)*rm1;
		}
	}

	is1=ia1*is1+ic1; is1=is1-m1*(is1/m1);
	is2=ia2*is2+ic2; is2=is2-m2*(is2/m2);
	is3=ia3*is3+ic3; is3=is3-m3*(is3/m3);

/*	Select a random entry from RAN and store a new number in its place */
	i=(ish*is3)/m3;
	r1=ran[i];
	ran[i]=(is1+is2*rm2)*rm1;
	return r1;
}
