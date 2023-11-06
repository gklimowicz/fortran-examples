/*	To calculate the Bessel function of first and second kind of order one
	for real positive argument
	For XB<=0 the function of second kind is not defined and function
	will return zero value  without any error message or flag.

	XB : (input) Argument at which the functions need to be evaluated
	BJ1 : (output) Calculated value of Bessel function of first kind
	BY1 : (output) Calculated value of Bessel function of second kind

	Required functions : None
*/
 
#include <math.h>

#define PI 3.14159265358979324
#define EUGAM 0.5772156649015328606    /*	the Euler's constant */

void bjy1(double xb, double *bj1, double *by1)

{
	double x,y,fn,fd,f,p1,q1,x2;
 
/*	Coefficients of rational function approximations */
	double a[7] = {1.156878510067059849e-02,  6.749406787503859073e-05,
                       2.614560837317535451e-07,  7.408815126464007290e-10,
                       1.577637796406197189e-12,  2.432945305413635549e-15,
                       2.257446839754248784e-18};
	double b[8] = {5.000000000000000074e-01, -5.671560744966475746e-02,
                       1.914864631812345532e-03, -2.821407888958585592e-05,
                       2.103168789640803591e-07, -8.322474383730280556e-10,
                       1.678871778708754849e-12, -1.372424374400306547e-15};

	double a0[7] = {1.186963690270342970e-02,  7.123839029323002736e-05,
                        2.848196759538669928e-07,  8.365303089083305885e-10,
                        1.857096246589275584e-12,  3.012506935004947491e-15,
                        2.996951174746838817e-18};
	double b0[8] = {2.500000000000000299e-01, -7.515759077432437273e-02,
                        3.430771992327672576e-03, -6.022315614557372919e-05,
                        5.067136874996839630e-07, -2.197514674456554803e-09,
                        4.768619679411702014e-12, -4.139491442515065355e-15};
    
	double a1[5] = { 8.659888261699365129e+01,
                         1.932665751369749084e+03,  1.172714583536277145e+04,
                         1.256737699073784218e+04, -6.147124347503755010e+02};
	double b1[5] = { 1.000000000000000011e+00,
                         8.671607011699346720e+01,  1.942669862370300601e+03,
                         1.194181952104744095e+04,  1.371467864585746530e+04};

	double a2[5] = { 1.021472573795463627e+02,  2.807865400111916226e+03,
                         2.280402060738415865e+04,  4.121116954504273053e+04,
                         3.501974669280301705e+03};
	double b2[5] = { 3.749999999999999461e-01,  3.820268245483084309e+01,
                         1.042753017477090289e+03,  8.289951986135169400e+03,
                         1.371889615877945967e+04};
 
	x=fabs(xb);

	if(x<8.0) {
/*	Use rational function approximation */
		y=x*x;
		fn=((((((b[7]*y+b[6])*y+b[5])*y+b[4])*y+b[3])*y+b[2])*y+b[1])*y+b[0];
		fd=((((((a[6]*y+a[5])*y+a[4])*y+a[3])*y+a[2])*y+a[1])*y+a[0])*y+1;
		*bj1=x*fn/fd;

		fn=((((((b0[7]*y+b0[6])*y+b0[5])*y+b0[4])*y+b0[3])*y+b0[2])*y+b0[1])*y+b0[0];
		fd=((((((a0[6]*y+a0[5])*y+a0[4])*y+a0[3])*y+a0[2])*y+a0[1])*y+a0[0])*y+1;
		if(x>0.0) *by1=2.*(*bj1*(log(x/2)+EUGAM)-1.0/x-x*fn/fd)/PI;
	}

	else {
/*	Use rational function approximations for P_1 and Q_1 */
		y=1./(x*x);
		fn=(((b1[4]*y+b1[3])*y+b1[2])*y+b1[1])*y+b1[0];
		fd=((((a1[4]*y+a1[3])*y+a1[2])*y+a1[1])*y+a1[0])*y+1;
		p1=fn/fd;

		fn=(((b2[4]*y+b2[3])*y+b2[2])*y+b2[1])*y+b2[0];
		fd=((((a2[4]*y+a2[3])*y+a2[2])*y+a2[1])*y+a2[0])*y+1;
		q1=fn/(fd*fabs(x));
 
		y=fabs(x);
		x2=y-0.75*PI;
		*by1=sqrt(2.0/(PI*y))*(p1*sin(x2)+q1*cos(x2));
		*bj1=sqrt(2.0/(PI*y))*(p1*cos(x2)-q1*sin(x2));
	}

	if(xb<=0.0) *by1=0.0;
	return;
}
