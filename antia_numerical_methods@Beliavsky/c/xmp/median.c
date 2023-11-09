/*	To calculate median, mean and higher moments of a distribution */

#include <stdio.h>
#include <math.h>

void shsort(double a[], int n);

main()
{
	int i,n;
	double amed,amean,sig,skew,aker,r,r2;
	double x[21]={49,73,46,45,67,47,28,49,61,71,49,57,52,79,40,41,55,69,41,71,40};

	n=21; shsort(x,n);
	amed=x[(n-1)/2];
	amean=0.0;
	for(i=0; i<n; ++i) amean=amean+x[i];
	amean=amean/n;
	sig=0.0;
	for(i=0; i<n; ++i) {
		r=(x[i]-amean); sig=sig+r*r;
	}
	sig=sqrt(sig/(n-1));
	skew=0.0; aker=0.0;
	for(i=0; i<n; ++i) {
		r=(x[i]-amean)/sig; r2=r*r;
		skew=skew+r2*r; aker=aker+r2*r2;
	}
	skew=skew/(n-1); aker=aker/(n-1)-3;
	printf(" Median = %f  Mean = %f  Standard deviation = %f \n",amed,amean,sig);
	printf(" Skewness = %f  Kurtosis = %f \n",skew,aker);
}

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

