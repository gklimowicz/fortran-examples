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
