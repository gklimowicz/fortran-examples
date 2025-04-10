/*	To calculate error function and complementary error function */

#include <stdio.h>
#include <math.h>

double erf(double xg);
double erfc(double xg);

main()
{
	int i,i1,j,k,n;
	double x,f1,f2;

	for(i1=0; i1<99; ++i1) {
		printf("type x,        (quits when x=-9999)\n");
		scanf(" %le", &x);
		if(x==-9999) return 0;

		f1=erf(x);
		f2=erfc(x);
		printf(" x= %e    erf(x) = %e     erfc(x) = %e\n",x,f1,f2);
	}
	return;
}



/*	To calculate the Error function for a real argument 

	Required functions : None
*/

#include <math.h>

double erf(double x)

{
	double fn,fd,f,y;

/*	The coefficients of rational function approximations */
	double a[7] = {4.891837297874833514e-01,  1.110175596090841322e-01,
                       1.526977188817169289e-02,  1.388143322498740953e-03,
                       8.446154421878829637e-05,  3.239915935631695227e-06,
                       6.200069065781009292e-08};
	double b[8] = {1.128379167095512575e+00,  1.758583405424390318e-01,
                       5.411290829613803886e-02,  3.805760982281134630e-03,
                       4.314620532020313106e-04,  1.423201530735308681e-05,
                       6.188698151904667110e-07,  2.293915472550064153e-09};
	double a1[7] = {3.040844316651787557e+01,  3.358927547920980388e+02,
                        1.703170048554673596e+03,  4.133195606956137156e+03,
                        4.555611776312694034e+03,  1.935778559646575488e+03,
                        2.051089518116697193e+02};
	double b1[8] = {5.641895835477562749e-01,  1.687403209467950089e+01,
                        1.813522721872712655e+02,  8.779664433851135986e+02,
                        1.965115658619443782e+03,  1.865558781108286245e+03,
                        5.828865035607128761e+02,  2.558559157228883880e+01};
 
	if(fabs(x)<2.0) {
/*	Use approximation for Erf(x)/x */
		y=x*x;
		fn=((((((b[7]*y+b[6])*y+b[5])*y+b[4])*y+b[3])*y+b[2])*y+b[1])*y+b[0];
		fd=((((((a[6]*y+a[5])*y+a[4])*y+a[3])*y+a[2])*y+a[1])*y+a[0])*y+1;
		f=x*fn/fd;
	}

	else {
/*	Use approximation for x*Exp(x^2)*Erfc(x) */
		y=1.0/(x*x);
		fn=((((((b1[7]*y+b1[6])*y+b1[5])*y+b1[4])*y+b1[3])*y+b1[2])*y+b1[1])*y+b1[0];
		fd=((((((a1[6]*y+a1[5])*y+a1[4])*y+a1[3])*y+a1[2])*y+a1[1])*y+a1[0])*y+1;
		f=1.-exp(-x*x)*fn/(fd*fabs(x));
		if(x<0.0) f=-f;
	}
	return f;
}


/*	To calculate the complementary Error function for a real argument 

	Required functions : None
*/
 
#include <math.h>

double erfc(double x)

{
	double fn,fd,f,y;

/*	The coefficients of rational function approximations */
	double a[7] = {4.891837297874833514e-01,  1.110175596090841322e-01,
                       1.526977188817169289e-02,  1.388143322498740953e-03,
                       8.446154421878829637e-05,  3.239915935631695227e-06,
                       6.200069065781009292e-08};
	double b[8] = {1.128379167095512575e+00,  1.758583405424390318e-01,
                       5.411290829613803886e-02,  3.805760982281134630e-03,
                       4.314620532020313106e-04,  1.423201530735308681e-05,
                       6.188698151904667110e-07,  2.293915472550064153e-09};
	double a1[7] = {3.040844316651787557e+01,  3.358927547920980388e+02,
                        1.703170048554673596e+03,  4.133195606956137156e+03,
                        4.555611776312694034e+03,  1.935778559646575488e+03,
                        2.051089518116697193e+02};
	double b1[8] = {5.641895835477562749e-01,  1.687403209467950089e+01,
                        1.813522721872712655e+02,  8.779664433851135986e+02,
                        1.965115658619443782e+03,  1.865558781108286245e+03,
                        5.828865035607128761e+02,  2.558559157228883880e+01};
 
	if(fabs(x)<2.0) {
/*	Use approximation for Erf(x)/x */
		y=x*x;
		fn=((((((b[7]*y+b[6])*y+b[5])*y+b[4])*y+b[3])*y+b[2])*y+b[1])*y+b[0];
		fd=((((((a[6]*y+a[5])*y+a[4])*y+a[3])*y+a[2])*y+a[1])*y+a[0])*y+1;
		f=1.-x*fn/fd;
	}

	else {
/*	Use approximation for x*Exp(x^2)*Erfc(x) */
		y=1.0/(x*x);
		fn=((((((b1[7]*y+b1[6])*y+b1[5])*y+b1[4])*y+b1[3])*y+b1[2])*y+b1[1])*y+b1[0];
		fd=((((((a1[6]*y+a1[5])*y+a1[4])*y+a1[3])*y+a1[2])*y+a1[1])*y+a1[0])*y+1;
		f=exp(-x*x)*fn/(fd*fabs(x));
		if(x<0.0) f=2.-f;
	}
	return f;
}
