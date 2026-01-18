/*	Real zero of a continuous function using bisection

	X : (output) Computed value of the zero using interpolation
	XL : (input/output) Lower limit of interval containing the zero
	XU : (input/output) Upper limit of interval containing the zero
		These limits will be replaced by refined limits after execution
	NB : (input) Number of bisections to be performed
	F : (input) Name of the function routine to calculate the function

	Error status is returned by the value of the function BISECT.
		0 value implies successful execution
		-1 implies that function value is zero at some trial point
		401 implies NB<=0 in which case no calculations are done
		421 implies that function has same sign at both limits
		and it is not possible to apply bisection

	Function F(X) must be supplied by the user.

	Required functions : F
*/

#include <math.h>

int bisect(double *x, double *xl, double *xu, int nb, double (*f) (double ))

{
	int k;
	double fl,fu,fx;

	if(nb<=0) return 401;

	fl=f(*xl); fu=f(*xu);
/*	If the function value is zero at either point then quit */
	if(fl==0.0) {*x=(*xl); return -1;}
	else if(fu==0) {*x=(*xu); return -1;}

/*	If the function has the same sign at both end points then quit */
	if((fl>0.0) == (fu>0.0)) return 421;

/*	Loop for bisection */
	for(k=1; k<=nb; ++k) {
		*x=(*xl+(*xu))/2.;
		fx=f(*x);
/*	If function is zero then quit */
		if(fx == 0.0) return -1;
		if((fx>0.0) == (fu>0.0)) {*xu=(*x); fu=fx;}
		else {*xl=(*x); fl=fx;}
	}

/*	linear interpolation between XL and XU */
	*x=((*xl)*fu-(*xu)*fl)/(fu-fl);
	return 0;
}
