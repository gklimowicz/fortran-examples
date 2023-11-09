/*	To generate uniformly distributed random numbers in interval (0,1)

	ISEED : (input/output) is an integer value used as the seed
		It should be initialised to negative value before first call
		and should not be modified between successive calls.

	Required functions : None
*/

#include <math.h>


double ranf(int *iseed)

{
	int m1=714025, ia1=1366, ic1=150889;
	int m2=214326, ia2=3613, ic2=45289;
	int m3=139968, ia3=3877, ic3=29573;
	int i,j, ish=43;
	double r1;

	static double rm1,rm2,ran[43];
	static int iflg=0, is1,is2,is3;
	
/*	Initialise on first call or when ISEED<0 */
	if(*iseed < 0 || iflg == 0) {
		iflg=1;
		rm1=1.0/m1;
		rm2=1.0/m2;

/*	Seeds for the three random number generators */
		is1=-(*iseed); is1=is1-m1*(is1/m1);
		is2=ia1*is1+ic1; is2=is2-m1*(is2/m1);
		is3=ia2*is2+ic2; is3=is3-m2*(is3/m2);
		*iseed=1;

/*	Store ISH random numbers in the array RAN */
		for(j=0; j<ish; ++j) {
			is1=ia1*is1+ic1; is1=is1-m1*(is1/m1);
			is2=ia2*is2+ic2; is2=is2-m2*(is2/m2);
			ran[j]=(is1+is2*rm2)*rm1;
		}
	}

	is1=ia1*is1+ic1; is1=is1-m1*(is1/m1);
	is2=ia2*is2+ic2; is2=is2-m2*(is2/m2);
	is3=ia3*is3+ic3; is3=is3-m3*(is3/m3);

/*	Select a random entry from RAN and store a new number in its place */
	i=(ish*is3)/m3;
	r1=ran[i];
	ran[i]=(is1+is2*rm2)*rm1;
	return r1;
}
