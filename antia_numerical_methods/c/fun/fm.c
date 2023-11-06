/*	To calculate the error in rational function approximation
	For use with function REMES. It is called by BRENTM to find
	extrema of error curve. It may not be used for any other purpose.

	X : (input) the value at which error is to be calculated
	FM = (FUN(X) - FUND(X)*R_mk(X))*SI  is the calculated error

	The parameters for rational function are passed through global variables
	MM,KK are degree of numerator and denominator
	AA is an array containing the coefficient of rational function approximation
	SI is -1 for maximum and +1 for minimum. The error is multiplied
		by SI to make it a minimum.
		For initial scan SI>10 and function is not evaluated.

	Functions FUN(x) and FUND(x) must be provided by the user to
		seek approximation of form FUN(X) = FUND(X)*RMK(X)

	Required functions : FUN, FUND
*/

/* To pass on these parameters to FM put these statements before that */
#define NMAXR 50
int MM, KK;
double AA[NMAXR],SI;

double fm(double x)

{
	int j,nk;
	double fd,fn,f;

	nk=MM+KK;
	f=0.0;
	if(SI<10.0) f=fun(x);

/*	Calculate the numerator using nested multiplication */
	fn=AA[KK+MM];
	for(j=1; j<=MM; ++j) fn=fn*x+AA[nk-j];

	if(KK>0) {
/*	Calculate the denominator using nested multiplication */
		fd=AA[KK-1];
		for(j=2; j<=KK; ++j) fd=fd*x+AA[KK-j];
		fd=fd*x+1.0;
	}
	else fd=1.0;

	f=f-fund(x)*fn/fd;
	if(SI<0.0) f=-f;
	return f;
}
