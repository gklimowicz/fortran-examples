/*	To find cascade sum of a series

	TERM : (input) Name of function to calculate the Ith term
	      It is possible to supply the terms in an array TERM
	      if TERM is declared as an array. (see CASSUM_A)
	N : (input) Number of terms to be summed

	Function TERM(I) to calculate the Ith term must be supplied by
	the user

	Required functions : TERM
*/

#include <math.h>


double cassum(double (*term)(int ), int n)

{
	int i, j, q[30], n2max=30;
	double s[30], sum;
	enum boolean {NO, YES};

	sum=0.0;
	for(i=0; i<n2max; ++i)  q[i]=YES;

	for(i=0; i<n; ++i) {
		if(q[0]) {s[0]=term(i); q[0]=NO;}
/*	If a pair is formed add the sum to higher level in the binary tree */
		else {
			s[0]=s[0]+term(i); q[0]=YES;
			for(j=1; j<n2max; ++j) {
				if(q[j]) {s[j]=s[j-1]; q[j]=NO; break;}
				else {s[j]=s[j]+s[j-1]; q[j]=YES;}
			}
			if(j>=n2max) sum=sum+s[n2max-1];
		}
	}

/*	Find the sum by adding up the incomplete pairs */

	for(j=0; j<n2max-1; ++j) {
		if(q[j]==NO) {
			if(q[j+1]) {s[j+1]=s[j]; q[j+1]=NO;}
			else {s[j+1]=s[j+1]+s[j];}
		}
	}
	if(q[n2max-1]==NO) sum=sum+s[n2max-1];
	return sum;
}
