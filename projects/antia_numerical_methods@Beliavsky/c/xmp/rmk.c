/*	To evaluate rational function and its derivative */

#include <stdio.h>
#include <math.h>

double rmk(int m, int k, double a[], double b[], double x);
double rmkd(int m, int k, double a[], double b[], double x, double *df);

main()
{
	int i,i1,j,k,n;
	double x,f1,f2,df1,df2,a[50],b[50];

	printf("type  n, k = degree of numerator and denominator \n");
	scanf(" %d %d",&n,&k);
	printf("type coefficients in numerator, starting with highest degree term\n");
	for(i=n; i>=0; --i) scanf(" %le",&a[i]);
	printf("  coefficients in numerator \n");
	for(i=n; i>=0; --i) printf(" %e",a[i]);
	printf("\n");

	printf("type coefficients in denominator, starting with highest degree term\n");
	for(i=k; i>=0; --i) scanf(" %le",&b[i]);
	printf("  coefficients in denominator \n");
	for(i=k; i>=0; --i) printf(" %e",b[i]);
	printf("\n");

	for(i1=0; i1<99; ++i1) {
		printf("type  x         (quits when x<-10000)\n");
		scanf(" %le", &x);
		if(x<-10000) return 0;

		f1=rmk(n,k,a,b,x);
		printf(" n = %d   k = %d    x= %e   rnk =  %e\n",n,k,x,f1);
		f2=rmkd(n,k,a,b,x,&df1);
		printf(" rnk = %e  rnk' = %e \n",f2,df1);
	}
	return;
}



/*	To evaluate a rational function at any point

	M : (input) Degree of numerator
	K : (input) Degree of denominator
	A : (input) Array of length M+1 containing the coefficients
		of the polynomial in numerator. A[0] is the constant term
		while A[M] is the coefficient of X**M.
	B : (input) Array of length K+1 containing the coefficients
		of the polynomial in denominator. B[0] is the constant term
		while B[K] is the coefficient of X**K.
	X : (input) The value of X at which the rational function needs
		to be evaluated
	The value of rational function will be returned as RMK

	Required functions : None
*/

double rmk(int m, int k, double a[], double b[], double x)

{
	int i;
	double den, rnu;

	rnu=a[m];
	for(i=m-1; i>=0; --i) rnu=rnu*x+a[i];

	den=b[k];
	for(i=k-1; i>=0; --i) den=den*x+b[i];

	return rnu/den;
}



/*	To evaluate a rational function and its derivative at any point


	M : (input) Degree of numerator
	K : (input) Degree of denominator
	A : (input) Array of length M+1 containing the coefficients
		of the polynomial in numerator. A[0] is the constant term
		while A[M] is the coefficient of X**M.
	B : (input) Array of length K+1 containing the coefficients
		of the polynomial in denominator. B[0] is the constant term
		while B[K] is the coefficient of X**K.
	X : (input) The value of X at which the rational function needs
		to be evaluated
	DF : (output) The first derivative of the rational function at X
	The value of rational function will be returned as RMKD

	Required functions : None
*/
 
double rmkd(int m, int k, double a[], double b[], double x, double *df)

{
	int i;
	double den, rnu, den1, dnu;

	rnu=a[m]; dnu=0.0;
	for(i=m-1; i>=0; --i) {dnu=dnu*x+rnu; rnu=rnu*x+a[i];}

	den=b[k]; den1=0.0;
	for(i=k-1; i>=0; --i) {den1=den1*x+den; den=den*x+b[i];}

	*df=dnu/den-rnu*den1/(den*den);
	return rnu/den;
}
