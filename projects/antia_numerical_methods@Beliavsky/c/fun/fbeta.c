/*	To calculate the integrand for incomplete beta function
	The variable AA, BB need to be defined outside

	Required routines : none
*/

double AA, BB;

#include <math.h>

	double fbeta(double x)
{
	return pow(x,AA-1)*pow(1-x,BB-1);
}

