/*	To balance a general real matrix by exact diagonal similarity transformation

	A : (input/output) Array of length IA*N containing the matrix
		elements. After execution it will contain the balanced matrix
	N : (input) Order of matrix
	IA : (input) Second dimension of array A as declared in the calling
		function
	B : (input) Base of the arithmetic to be used for balancing.
		Should generally be 2.
	LOW : (output) Index of column such that in balanced matrix A[i][j]=0
		if i>j and j<LOW, the first LOW eigenvalues are isolated
	IGH : (output) Index of row such that in balanced matrix A[i][j]=0
		if i>j and i>IGH, the last N-IGH-1 eigenvalues are isolated
		After balancing only the sub-matrix between rows and columns
		LOW to IGH need to be considered
	D : (output) Array of length N containing information about
		transformations used for balancing. D[LOW] to D[IGH]
		contain elements used for balancing, while other elements
		contain information about permutations used

	Error status is returned by the value of the function BALANC.
		0 value implies successful execution
		112 implies that N<=1 or N>IA

	Required functions : None
*/

#include <math.h>

int balanc(double *a, int n, int ia, double b, int *low, int *igh, double d[])

{
	int i,j,qc;
	double b2,r,t,c,f,g,s, gam=0.95;

	if(n<=1 || n>ia) return 112;
	b2=b*b;
	*low=0; *igh=n-1;

/*	Search for rows isolating an eigenvalue */
	do {
		for(j=(*igh); j>=0; --j) {
			r=0.0;
			for(i=0; i<=j-1; ++i) r=r+fabs(a[i+j*ia]);
			for(i=j+1; i<=(*igh); ++i) r=r+fabs(a[i+j*ia]);

			if(r==0.0) {
/*	Push the row to bottom */
				d[*igh]=j;
				if(j != (*igh)) {
					for(i=0; i<=(*igh); ++i) {
						t=a[j+i*ia];
						a[j+i*ia]=a[*igh+i*ia];
						a[*igh+i*ia]=t;
					}
					for(i=(*low); i<n; ++i) {
						t=a[i+j*ia];
						a[i+j*ia]=a[i+(*igh)*ia];
						a[i+(*igh)*ia]=t;
					}
				}
				*igh=(*igh)-1;
				break;
			}
		}
	} while(r==0.0);
			
/*	Search for columns isolating an eigenvalue */
	do {
		for(j=(*low); j<=(*igh); ++j) {
			c=0.0;
			for(i=(*low); i<=j-1; ++i) c=c+fabs(a[j+i*ia]);
			for(i=j+1; i<=(*igh); ++i) c=c+fabs(a[j+i*ia]);

/*	Move the column to the left end */
			if(c==0) {
				d[*low]=j;
				if(j != (*low)) {
					for(i=0; i<=(*igh); ++i) {
						t=a[j+i*ia];
						a[j+i*ia]=a[*low+i*ia];
						a[*low+i*ia]=t;
					}
					for(i=(*low); i<n; ++i) {
						t=a[i+j*ia];
						a[i+j*ia]=a[i+(*low)*ia];
						a[i+(*low)*ia]=t;
					}
				}
				*low=(*low)+1;
				break;
			}
		}
	} while(c==0);
				
/*	Balance the submatrix in rows LOW to IGH */
	for(i=(*low); i<=(*igh); ++i) d[i]=1;

	do {
		qc=0;
		for(i=(*low); i<=(*igh); ++i) {
			c=0.0; r=0.0;
			for(j=(*low); j<=(*igh); ++j) {
				if(j != i) {
					c=c+fabs(a[i+j*ia]);
					r=r+fabs(a[j+i*ia]);
				}
			}
			g=r/b; f=1.0; s=c+r;

			while (c<g) {f=f*b; c=c*b2;}
			g=r*b;

			while(c>=g) {f=f/b; c=c/b2;}

/*	Apply the transformation */
			if((c+r)/f < gam*s) {
				g=1./f;
				d[i]=d[i]*f;
				qc=1;
				for(j=(*low); j<n; ++j) a[j+i*ia]=a[j+i*ia]*g;
				for(j=0; j<=(*igh); ++j) a[i+j*ia]=a[i+j*ia]*f;
			}
		}
	} while(qc==1);
	return 0;
}
