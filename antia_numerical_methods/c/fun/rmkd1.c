/*	To evaluate a rational function and its derivative at any point,
	the constant term in denominator is assumed to be unity.

	M : (input) Degree of numerator
	K : (input) Degree of denominator
	A : (input) Array of length M+1 containing the coefficients
		of the polynomial in numerator. A[0] is the constant term
		while A[M] is the coefficient of X**M.
	B : (input) Array of length K containing the coefficients
		of the polynomial in denominator. B[J-1] is the coefficient of X**J.
		The constant term is assumed to be unity and is not supplied.
	X : (input) The value of X at which the rational function needs
		to be evaluated
	DF : (output) The first derivative of the rational function at X
	The value of rational function will be returned as RMKD1

	Required functions : None
*/
 
double rmkd1(int m, int k, double a[], double b[], double x, double *df)

{
	int i;
	double den, rnu, den1, dnu;

	rnu=a[m]; dnu=0.0;
	for(i=m-1; i>=0; --i) {dnu=dnu*x+rnu; rnu=rnu*x+a[i];}

	den=b[k-1]; den1=0.0;
	for(i=k-2; i>=0; --i) {den1=den1*x+den; den=den*x+b[i];}
	den1=den1*x+den; den=den*x+1.0;

	*df=dnu/den-rnu*den1/(den*den);
	return rnu/den;
}