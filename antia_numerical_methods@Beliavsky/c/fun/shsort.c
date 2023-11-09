/*       To sort an array of length n in ascending order, using shell sort algorithm

       a : (input) Real array of length n, which is to be sorted
       n : (input) No. of elements of array a to be sorted

       Required routines : None
*/

#include <math.h>

void shsort(double a[], int n)
{
	int i,j,i1,l,m;
	double an,al,t;

	an=n; l=log(an)/log(2.0)+1; al=l;
	m=pow(2.0,al);
	if(m>=2*n) m=m/2;
	while(m>0)
	{
		m=m/2;
		for(i=m; i<n; i++)
		{
			t=a[i]; i1=i;
			while(a[i1-m]>t && i1>=m)
			{
				a[i1]=a[i1-m]; i1=i1-m;
			}
			a[i1]=t;
		}
	}
}
