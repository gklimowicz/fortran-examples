/*	To generate uniformly distributed random numbers in interval (0,1)

	SEED : (input/output) is a real value used as the seed
		It should be positive during initial call and
		should not be modified between different calls

	Required functions : None
*/

#include <math.h>

double ran1(double *seed)

{
	double am=2147483648e0, a=45875e0, ac=453816693e0, an=2147483647e0;

	*seed=fmod((*seed)*a+ac,am);
	return (*seed)/an;
}
	
