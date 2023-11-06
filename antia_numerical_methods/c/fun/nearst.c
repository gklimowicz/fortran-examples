/*	To locate the nearest point in an ordered table using bisection

	XB : (input) given value of x for which nearest point is needed
	X : (input) array of length NTAB containing table of values
	NTAB : (input) length of table
	After execution X[NEARST] is the tabular point closest to XB 

	Required functions : None
*/

#include <math.h>

int nearst(double xb, double x[], int ntab)

{
	int low, igh, mid;

	low=0; igh=ntab-1;
	if((xb < x[low]) != (xb < x[igh]) ) {

/*	If the point is within the range of table, then locate it by bisection */

		while(igh-low > 1) {
			mid=(low+igh)/2;
			if((xb < x[mid]) == (xb < x[low])) low=mid;
			else igh=mid;
		}
	}

	if(fabs(xb-x[low]) < fabs(xb-x[igh])) return low;
	else return igh;
}
				
