/*	To calculate the Bessel function of first and second kind of order zero
	for real positive argument
	For XB<=0 the function of second kind is not defined and function
	will return zero value  without any error message or flag

	XB : (input) Argument at which the functions need to be evaluated
	BJ0 : (output) Calculated value of Bessel function of first kind
	BY0 : (output) Calculated value of Bessel function of second kind

	Required functions : None
*/

#include <math.h>

#define PI 3.14159265358979324
#define EUGAM 0.5772156649015328606    /*	the Euler's constant */

void bjy0(double xb, double *bj0, double *by0)

{
	double x,y,fn,fd,f,p0,q0;
 
/*	Coefficients of rational function approximations */
	double a[8] = {1.293686560051304152e-02,  8.573459862295151747e-05,
                       3.854769244046149702e-07,  1.308534328117880493e-09,
                       3.512360907188715842e-12,  7.512575042421009221e-15,
                       1.229302278444845702e-17,  1.311883486088925264e-20};
	double b[8] = {9.999999999999999878e-01, -2.370631343994868513e-01,
                       1.247651819849453565e-02, -2.529374255010058573e-04,
                       2.411267406461247155e-06, -1.159484705672466498e-08,
                       2.730546745501229851e-11, -2.517936655103065990e-14};

	double a0[7] = {1.089079731266387424e-02,  5.954632605213292419e-05,
                        2.150109922530480401e-07,  5.641082188778387960e-10,
                        1.102341761343675716e-12,  1.539990321465010920e-15,
                        1.263081729204829936e-18};
	double b0[8] = {2.500000000000000006e-01, -2.071480067183403591e-02,
                        5.553511120900719150e-04, -6.804373640943337406e-06,
                        4.346149688717712144e-08, -1.505635199331021665e-10,
                        2.703193275976574669e-13, -1.993047807317608951e-16};

	double a1[5] = { 8.911849018950665793e+01,
                         2.078818787053760203e+03,  1.366258799766718466e+04,
                         1.800383785973922830e+04,  4.923440494847201509e+02};
	double b1[5] = { 9.999999999999999908e-01,
                         8.904817768950681616e+01,  2.072664795311476688e+03,
                         1.352584337655551999e+04,  1.723138433448795889e+04};

	double a2[5] = { 1.046366195300779895e+02,  2.980727727381642723e+03,
                         2.570418404044668245e+04,  5.291014161889741749e+04,
                         9.497228391199055149e+03};
	double b2[5] = {-1.249999999999999534e-01, -1.300633525376058087e+01,
                        -3.651542590150084093e+02, -3.016744074771875522e+03,
                        -5.251679479249748063e+03};
 
	x=fabs(xb);

	if(x<8.0) {
/*	Use rational function approximation */
		y=x*x;
		fn=((((((b[7]*y+b[6])*y+b[5])*y+b[4])*y+b[3])*y+b[2])*y+b[1])*y+b[0];
		fd=(((((((a[7]*y+a[6])*y+a[5])*y+a[4])*y+a[3])*y+a[2])*y+a[1])*y+a[0])*y+1;
		*bj0=fn/fd;

		fn=((((((b0[7]*y+b0[6])*y+b0[5])*y+b0[4])*y+b0[3])*y+b0[2])*y+b0[1])*y+b0[0];
		fd=((((((a0[6]*y+a0[5])*y+a0[4])*y+a0[3])*y+a0[2])*y+a0[1])*y+a0[0])*y+1;
		if(x>0.0) *by0=2.*(*bj0*(log(x/2)+EUGAM)+y*fn/fd)/PI;
	}

	else {
/*	Use rational function approximations for P_0 and Q_0 */
		y=1./(x*x);
		fn=(((b1[4]*y+b1[3])*y+b1[2])*y+b1[1])*y+b1[0];
		fd=((((a1[4]*y+a1[3])*y+a1[2])*y+a1[1])*y+a1[0])*y+1;
		p0=fn/fd;

		fn=(((b2[4]*y+b2[3])*y+b2[2])*y+b2[1])*y+b2[0];
		fd=((((a2[4]*y+a2[3])*y+a2[2])*y+a2[1])*y+a2[0])*y+1;
		q0=fn/(fd*fabs(x));
 
		y=fabs(x);
		*by0=sqrt(2.0/(PI*y))*(p0*sin(y-PI/4.0)+q0*cos(y-PI/4.0));
		*bj0=sqrt(2.0/(PI*y))*(p0*cos(y-PI/4.0)-q0*sin(y-PI/4.0)); 
	}

	if(xb<=0.0) *by0=0.0;
	return;
}
