/*	To locate a given point between two points of an ordered table

	XB : (input) The point which needs to be located
	X : (input) Array of length NP containing the ordered table
	NP : (input) Length of table
	LOCATE should give the value such that XB is between X[LOCATE] and 
		X[LOCATE+1]

	Required functions : None
*/

#include <math.h>

int locate(double xb, double x[], int np)

{
	int low, igh, mid; 

	low=0;
	igh=np-1;
	if( (xb<x[low]) == (xb<x[igh]) ) {
		if(fabs(xb-x[low]) > fabs(xb-x[igh])) low=igh-1;
	}
	else {

/*	If the point is within the range of table locate it by bisection */
		while(igh-low > 1) {
			mid=(low+igh)/2;
			if( (xb<x[mid]) == (xb<x[low])) low=mid;
			else igh=mid;
		}
	}
	return low;
}
