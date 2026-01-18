/*	To find eigenvalues and left-eigenvectors of general real matrix 
	Since real arithmetic is used only eigenvectors corresponding to real
	eigenvalues can be calculated correctly, however all eigenvalues are
	calculated including complex ones */

#include <stdio.h>
#include <math.h>

int balanc(double *a, int n, int ia, double b, int *low, int *igh, double d[]);
int balbak_l(int n, int low, int igh, double *cz, int m, int iz, double d[]);
int elmhes(double *a, int n, int ia, int low, int igh, int inc[]);
int hqr(double *h, int nn, int ih, double er[], double ei[], double reps);
int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);
int invit_l(double *a, int m, int ia, double *p, double u[], int iflg,
	double *ei, double *erc, double reps, int *nit);

main()
{
	int i,i1,j,n,low,igh, id, iflg, ier,np,inc[20];
	double hh, d[100],e[100], c[4][4],b,ei[20],eps,reps,p0,v[20],ev[4][4];
	double pi,pr;
	double a[4][4]={-2,2,2,2, -3,3,2,2, -2,0,4,2, -1,0,0,5};
/*	double a[7][7]={6,0,0,0,0,1,0, 0,4,0,3.0e-8,1.e-4,2.e-4,1.e-3,
		1,10000,7,0,0,-2,20, 0,2.e8,0,1,-40000,30000,-4.e5,
		-2,-30000,0,1.e-4,2,2,40, 0,0,0,0,0,0,0,
		0,1000,0,4.0e-5,0.1,-0.2,3};  */

	id=4; n=4; np=0; reps=1.e-300;
	b=2.0;
	printf(" The matrix is :\n");
	for(i=0; i<n; ++i) {
		for(j=0; j<n; ++j) printf(" %e ",a[i][j]);
		printf("\n");
	}

/*	Balance the matrix */
	i=balanc(&a[0][0],n,id,b,&low,&igh,d);
	printf(" ier = %d    n =  %d    low = %d   igh = %d \n",i,n,low,igh);
	printf(" diagonal transformation matrix for balancing :\n");
	for(i=0; i<n; ++i) printf(" %e ",d[i]);
	printf(" \n  Balanced matrix :\n");
	for(i=0; i<n; ++i) {
/*	Preserve the matrix for back transformation of eigenvectors */
		for(j=0; j<n; ++j) c[i][j]=a[i][j];
		for(j=0; j<n; ++j) printf(" %e ",a[i][j]);
		printf("\n");
	}

/*	Reduce the matrix to Hessenberg form */
	i=elmhes(&a[0][0],n,id,low,igh,inc);
	printf(" i= %d  interchanges = ",i);
	for(i=0; i<n; ++i) printf(" %d ",inc[i]);
	printf(" \n transformed Hessenberg matrix  :\n ");
	for(i=0; i<n; ++i) {
		for(j=0; j<n; ++j) printf(" %e ",a[i][j]);
		printf(" \n");
	}

/*	Find eigenvalues of Hessenberg matrix */
	reps=1.e-14;
	i=hqr(&a[0][0],n,id,e,ei,reps);
	printf(" ier = %d    eigenvalues : ",i);
	for(i=0; i<n; ++i)  printf("   %e %e \n",e[i],ei[i]);
	printf(" \n");

/*	Finding eigenvectors using inverse iteration */
	for(i=0; i<n; ++i) {
/*	Perturb the eigenvalue slightly to avoid zero pivot 
	The imaginary part is ignored so only eigenvectors for real eigenvalues 
	can be correctly calculated, */
		p0=e[i]*(1.+20.0*reps)+reps*reps;
		iflg=0;
		np=100;
		for(j=0; j<n; ++j) v[j]=1.0;
		ier=invit_l(&c[0][0],n,id,&p0,v,iflg,&pi,&pr,reps,&np);
		for(j=0; j<n; ++j) {
			ev[j][i]=v[j];
		}
	} 

/*	Back transform the eigenvectors of balanced matrix */
		i=balbak_l(n,low,igh,&ev[0][0],n,id,d);
		printf("   eigenvalues      eigenvectors \n");
		for(i=0; i<n; ++i) {
			printf(" %e    ",e[i]);
			hh=0.0;
			for(j=0; j<n; ++j) {
				if(fabs(ev[j][i])>hh) hh=fabs(ev[j][i]);
			}
			for(j=0; j<n; ++j) printf(" %e ",ev[j][i]/hh);
			printf(" \n");
		}
		

	return;
}



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




/*	Perform back-transformation of a set of left eigenvectors from
	those of balanced matrix to that for original matrix
	Only real eigenvectors are handled by this function
	Since transformation is real and linear, the real and imaginary
	parts of complex eigenvectors can be transformed by two separate calls
	to this function.

	N : (input) Order of matrix
	LOW, IGH : (input) After balancing only rows and columns from
		LOW to IGH need to be considered, since other rows or
		columns contain isolated eigenvalues
	CZ : (input/output) Array of length IZ*M containing the
		eigenvectors of balanced matrix. After execution it will
		be overwritten by eigenvectors of the original matrix.
	M : (input) Number of eigenvectors to be balanced
	IZ : (input) Second dimension of array CZ as declared in the calling
		function
	D : (input) Array of length N, containing information about
		transformation used for balancing.

	Returned value is always 0

	Required functions : None
*/

#include <math.h>

int balbak_l(int n, int low, int igh, double *cz, int m, int iz, double d[])

{
	int i,j,k;
	double s,cs;

	for(i=low; i<=igh; ++i) {
       s=1.0/d[i];  

	   for(j=0; j<m; ++j) cz[j+i*iz]=cz[j+i*iz]*s;
	}

	for(i=low-1; i>=0; --i) {
		k=d[i];
		if(k != i) {
/*	Exchange the corresponding rows */
			for(j=0; j<m; ++j) {
				cs=cz[j+i*iz]; cz[j+i*iz]=cz[j+k*iz]; cz[j+k*iz]=cs;
			}
		}
	}

	for(i=igh+1; i<n; ++i) {
		k=d[i];
		if(k != i) {
/*	Exchange the corresponding rows */
			for(j=0; j<m; ++j) {
				cs=cz[j+i*iz]; cz[j+i*iz]=cz[j+k*iz]; cz[j+k*iz]=cs;
			}
		}
	}
	return 0;
}



/*	To reduce a general real matrix to Hessenberg form using stabilised
	elementary transformations
	It is advisable to balance the matrix before applying this transformations

	A : (input/output) Array of length IA*N containing the matrix
		After execution the reduced matrix will be overwritten on
		the same array
	N : (input) Order of matrix
	IA : (input) Second dimension of array A as declared in the calling
		function
	LOW, IGH : (input) After balancing only rows and columns from
		LOW to IGH need to be considered, since other rows or
		columns contain isolated eigenvalues. If the matrix is
		not balanced use LOW=0, IGH=N-1
	INC : (output) Integer array of length N, which will contain information
		about row and column interchanges during reduction
		
	Error status is returned by the value of the function ELMHES.
		0 value implies successful execution
		113 implies that N<=1 or N>IA

	Required functions : None
*/

#include <math.h>

int elmhes(double *a, int n, int ia, int low, int igh, int inc[])

{
	int i,j,m;
	double amax,t;

	if(n<=1 || n>ia) return 113;

	if(low>igh-2) return 0;

	for(i=0; i<n; ++i) inc[i]=i;
	for(m=low+1; m<=igh-1; ++m) {
		i=m;
/*	Find the pivot */
		amax=0.0;
		for(j=m; j<=igh; ++j) {
			if(fabs(a[m-1+j*ia]) > fabs(amax)) {amax=a[m-1+j*ia]; i=j;}
		}
		inc[m]=i;

		if(i != m) {
/*	Interchange the corresponding rows and columns */
			for(j=m-1; j<n; ++j) {
				t=a[j+i*ia]; a[j+i*ia]=a[j+m*ia]; a[j+m*ia]=t;
			}
			for(j=0; j<=igh; ++j) {
				t=a[i+j*ia]; a[i+j*ia]=a[m+j*ia]; a[m+j*ia]=t;
			}
		}

		if(amax != 0.0) {
/*	Perform Gaussian elimination */
			for(i=m+1; i<=igh; ++i) {
				t=a[m-1+i*ia];
				if(t != 0.0) {
					t=t/amax;
					a[m-1+i*ia]=t;
					for(j=m; j<n; ++j) a[j+i*ia]=a[j+i*ia]-t*a[j+m*ia];
					for(j=0; j<=igh; ++j) a[m+j*ia]=a[m+j*ia]+t*a[i+j*ia];
				}
			}
		}
	}
	return 0;
}



/*	Solution of a system of linear equations using Gaussian elimination
	with partial pivoting

	N : (input) Number of equations to be solved
	NUM : (input) Number of different sets (each with N equations) of
	         equations to be solved
	A : (input/output) The matrix of coefficient of size LJ*N
	        A[i][j] is the coefficient of x_j in ith equation
	     	at output it will contain the triangular decomposition
	X : (input/output) The matrix containing right hand sides (size LJ*NUM)
	        X[j][i] is the ith element of jth right hand side
	     	at output it will contain the solutions
	DET : (output) The determinant of the matrix
	INC : (output) Integer array of length N containing information about
		interchanges performed during elimination
	LJ : (input) Second dimension of arrays A and X in calling function
	IFLG : (input) Integer parameter to specify the type of computation required
		If IFLG<=0, both elimination and solution are
			done and IFLG is set to 2
		If IFLG=1, only elimination is done and IFLG is set to 2
		If IFLG>=2 only solution is calculated, the triangular
		    decomposition should have been calculated earlier
		
	Error status is returned by the value of the function GAUELM.
		0 value implies successful execution
		101 implies (N<=0 or N>LJ) 
		121 implies some pivot turned out to be zero and hence
			matrix must be nearly singular

	Required functions : None
*/

#include <math.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg)

{
int i,j,k,km,l;
double r1,t1;

	if(n<=0 || n>lj) return 101;
 
	if((*iflg)<2) {
/*	Perform elimination  */
		*det=1.0;
		for(k=0; k<n-1; ++k) {
/*	Find the maximum element in the Kth column  */
			r1=0.0; km=k;
			for(l=k; l<n; ++l)
				if(fabs(a[l*lj+k])>r1) {r1=fabs(a[l*lj+k]); km=l;}

			inc[k]=km;
			if(km != k) {
/*	Interchange the rows if needed  */
				for(l=k; l<n; ++l) 
				{t1=a[k*lj+l]; a[k*lj+l]=a[km*lj+l]; a[km*lj+l]=t1;}
				*det=-(*det);
			}

			*det=(*det)*a[k*lj+k];
			if(a[k*lj+k]==0) return 121;

			for(l=k+1; l<n; ++l) {
				a[l*lj+k]=a[l*lj+k]/a[k*lj+k];
				for(i=k+1; i<n; ++i) a[l*lj+i]=a[l*lj+i]-a[l*lj+k]*a[k*lj+i];
			}
		}
		*det=(*det)*a[(n-1)*lj+n-1];
		inc[n-1]=n-1;
		if(a[(n-1)*lj+n-1]==0) return 121;

		if((*iflg)==1) {*iflg=2; return 0;}
		*iflg=2;
	}

/*	Solution for the num different right-hand sides  */
	for(j=0; j<num; ++j) {
/*	forward-substitution  */
		for(k=0; k<n-1; ++k) {
			if(k != inc[k])
			{t1=x[j*lj+k]; x[j*lj+k]=x[j*lj+inc[k]]; x[j*lj+inc[k]]=t1;}
			for(l=k+1; l<n; ++l) x[j*lj+l]=x[j*lj+l]-a[l*lj+k]*x[j*lj+k];
		}

/*	back-substitution  */

		x[j*lj+n-1]=x[j*lj+n-1]/a[(n-1)*lj+n-1];
		for(k=n-2; k>=0; --k) {
			for(l=n-1; l>=k+1; --l) x[j*lj+k]=x[j*lj+k]-x[j*lj+l]*a[k*lj+l];
			x[j*lj+k]=x[j*lj+k]/a[k*lj+k];
		}
	}
	return 0;
}



/*	To find all eigenvalues of an upper-Hessenberg matrix using QR algorithm

	H : (input/output) Array of length IH*NN containing the matrix
		elements. The contents of array are destroyed during execution
	NN : (input) Order of matrix
	IH : (input) Second dimension of array H as declared in the calling
		function
	ER : (output) Array of length NN containing the real part of
		the eigenvalues.
	EI : (output) Array of length NN containing the imaginary part of
		the eigenvalues.
	REPS : (input) Required relative accuracy in eigenvalues. It should
		be of the order of machine accuracy
		
	Error status is returned by the value of the function HQR.
		0 value implies successful execution
		114 implies that N<=1 or N>IH, in which case no
			calculations are done
		145 implies that QR iteration failed to converge at
			some stage and calculations are abandoned.

	Required functions : None
*/

#include <math.h>

int hqr(double *h, int nn, int ih, double er[], double ei[], double reps)

{
	int i,j,k,l,m,n,it, nit=30;
	double t,x,w,y,z,p,d,r,s,u;

	t=0.0; n=nn-1;
	if(nn<=1 || nn>ih) return 114;

/*	If all eigenvalues are found, then quit */
nexti: if(n<0) return 0;

/*	Loop for double QR iteration */
	for(it=1; it<=nit; ++it) {
/*	Search for a small subdiagonal element */
		for(l=n; l>=1; --l) {
			if(fabs(h[l-1+l*ih]) <= reps*(fabs(h[l-1+(l-1)*ih])+fabs(h[l+l*ih]))) break;
		}
		x=h[n+n*ih];

		if(l==n) {
/*	One eigenvalue is isolated */
			er[n]=x+t; ei[n]=0.0;
			n=n-1;
			goto nexti;
		}

		y=h[n-1+(n-1)*ih];
		w=h[n-1+n*ih]*h[n+(n-1)*ih];
		if(l == n-1) {
/*	A pair of eigenvalues are isolated */
			p=(y-x)/2.0;
			d=p*p+w;
			y=sqrt(fabs(d));
			x=x+t;
			if(d>0.0) {
/*	Pair of real eigenvalues */
				if(p<0.0) y=-y;
				y=p+y;
				er[n-1]=x+y; er[n]=x-w/y;
				ei[n-1]=0.0; ei[n]=0.0;
			}
			else {
/*	Pair of complex conjugate eigenvalues */
				er[n]=x+p; er[n-1]=er[n];
				ei[n-1]=y; ei[n]=-y;
			}
			n=n-2;
			goto nexti;
		}

		if(it==10 || it==20) {
/*	Apply special shifts */
			t=t+x;
			for(i=0; i<=n; ++i) h[i+i*ih]=h[i+i*ih]-x;
			s=fabs(h[n-1+n*ih])+fabs(h[n-2+(n-1)*ih]);
			x=0.75*s;
			y=x;
			w=-0.4375*s*s;
		}

/*	Search for two consecutive small subdiagonal elements */
		for(m=n-2; m>=l; --m) {
			z=h[m+m*ih];
			r=x-z; s=y-z;
			p=(r*s-w)/h[m+(m+1)*ih]+h[m+1+m*ih];
			u=h[m+1+(m+1)*ih]-z-r-s;
			r=h[m+1+(m+2)*ih];
			s=fabs(p)+fabs(u)+fabs(r);
			p=p/s; u=u/s; r=r/s;
			if(m==l) break;
			if(fabs(h[m-1+m*ih])*(fabs(u)+fabs(r))<= reps*fabs(p)*
				(fabs(h[m-1+(m-1)*ih])+fabs(z)+fabs(h[m+1+(m+1)*ih]))) break;
		}
		for(i=m+3; i<=n; ++i) {h[i-3+i*ih]=0.0; h[i-2+i*ih]=0.0;}
		h[m+(m+2)*ih]=0.0;

/*	Double QR transform for rows L to N and columns M to N */
		for(k=m; k<=n-1; ++k) {
			if(k != m) {
				p=h[k-1+k*ih];
				u=h[k-1+(k+1)*ih];
				r=0.0;
				if(k != n-1) r=h[k-1+(k+2)*ih];
				x=fabs(p)+fabs(u)+fabs(r);
				if(x != 0.0) {p=p/x; u=u/x; r=r/x;}
			}

			if(x != 0.0) {
				s=sqrt(p*p+u*u+r*r);
				if(p<0.0) s=-s;
				if(k != m) h[k-1+k*ih]=-s*x;
				else if(l != m) h[k-1+k*ih]=-h[k-1+k*ih];
				p=p+s; x=p/s; y=u/s; z=r/s; u=u/p; r=r/p;

/*	Row modification */
				for(j=k; j<=n; ++j) {
					p=h[j+k*ih]+u*h[j+(k+1)*ih];
					if(k != n-1) {
						p=p+r*h[j+(k+2)*ih];
						h[j+(k+2)*ih]=h[j+(k+2)*ih]-p*z;
					}
					h[j+(k+1)*ih]=h[j+(k+1)*ih]-p*y;
					h[j+k*ih]=h[j+k*ih]-p*x;
				}

/*	Column modification */
				j=k+3; if(n<j) j=n;
				for(i=l; i<=j; ++i) {
					p=x*h[k+i*ih]+y*h[k+1+i*ih];
					if(k != n-1) {
						p=p+z*h[k+2+i*ih];
						h[k+2+i*ih]=h[k+2+i*ih]-p*r;
					}
					h[k+1+i*ih]=h[k+1+i*ih]-p*u;
					h[k+i*ih]=h[k+i*ih]-p;
				}
			}
		}
	}

	return 145;
}


/*	Real eigenvalues and left eigenvector of a real matrix using inverse iteration

	A : (input) Array of length IA*M containing the matrix elements
	M : (input) Order of the matrix
	IA : (input) The second dimension of A as declared in the calling function, IA>=M
	P : (input/output) Initial value of the shift. This will be modified
		by the function if IFLG>0
	U : (input/output) Array of length M, which should specify the
		initial approximation to eigenvector. After execution it
		will contain the calculated left eigenvector.
	IFLG : (input) Integer variable to specify the type of iteration required
		If IFLG=0 the shift P is kept fixed
		If IFLG=1 the shift P is varied using Rayleigh quotient
		If IFLG=2 the shift P is varied using max(V_s+1)
	EI : (output) Estimated eigenvalue using simple inverse iteration
	ERC : (output) Estimated eigenvalue using Rayleigh quotient
	REPS : (input) Required absolute accuracy. Iteration is terminated
		when all components of eigenvector and the eigenvalue have
		converged to REPS.
	NIT : (input/output) Number of iterations required. If it is
		zero or negative NIT is set to NIT0=100
		
	Error status is returned by the value of the function INVIT_L.
		0 value implies successful execution
		106 implies that M<=1 or M>IA, in which case no
			calculations are done
		141 implies that vector is zero at some stage and
			calculations are aborted
		142 implies that inverse iteration has failed to converge

	Required functions : GAUELM
*/

#include <math.h>
#include <stdlib.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);

int invit_l(double *a, int m, int ia, double *p, double u[], int iflg,
	double *ei, double *erc, double reps, int *nit)

{
	int i,j,k,ifl,km,ier, nit0=100;
	double epi,r1,ukm,s1,s2,det;
	int *iwk;
	double *wk,*v;

	if(m<=1 || m> ia) return 106;

	iwk=(int *) calloc((size_t) m, sizeof(int));
	v=(double *) calloc((size_t) m, sizeof(double));
	wk=(double *) calloc((size_t) (m*m), sizeof(double));

/*	Copy the matrix to WK and apply the shift */
	for(i=0; i<m; ++i) {
		v[i]=u[i];
		for(j=0; j<m; ++j) wk[j+i*m]=a[i+j*m];
		wk[i+i*m]=a[i+i*m]-(*p);
	}

	ifl=1;
	if(*nit<=0) *nit=nit0;
/*	Perform Gaussian elimination on A-pI */
	ier=gauelm(m,1,wk,u,&det,iwk,m,&ifl);
	if(ier>0) {free(wk); free(v); free(iwk); return ier;}

	epi=0.0;
/*	Loop for inverse iteration */
	for(j=1; j<=(*nit); ++j) {
		ier=gauelm(m,1,wk,u,&det,iwk,m,&ifl);
		if(ier>0) {free(wk); free(v); free(iwk); return ier;}

/*	Normalising the vector U */
		r1=0.0; km=0;
		for(k=0; k<m; ++k) {
			if(r1<fabs(u[k])) {r1=fabs(u[k]); km=k;}
		}
		ukm=u[km];
/*	If the vector is zero, then quit */
		if(ukm==0.0) {free(wk); free(v); free(iwk); return 141;}

/*	The eigenvalue */
		*ei=v[km]/ukm+(*p);
		s1=0.0; s2=0.0;
/*	Calculating the Rayleigh quotient */
		for(k=0; k<m; ++k) {
			s1=s1+u[k]*v[k];
			s2=s2+fabs(u[k])*fabs(u[k]);
			u[k]=u[k]/ukm;
		}
		*erc=(*p)+s1/s2;

/*	Convergence check */
		r1=fabs(*ei-epi);
		for(i=0; i<m; ++i) {
			if(r1<fabs(v[i]-u[i])) r1=fabs(v[i]-u[i]);
			v[i]=u[i];
		}
		if(fabs(r1)<reps) {free(wk); free(v); free(iwk); return 0;}
		epi=(*ei);

		if(iflg>=1) {
/*	Update the shift */
			*p=(*erc);
			if(iflg==2) *p=(*ei);
/*	Setting up the new matrix A-pI */
			for(i=0; i<m; ++i) {
				for(k=0; k<m; ++k) wk[k+i*m]=a[i+k*m];
				wk[i+i*m]=a[i+i*m]-(*p);
			}
			ifl=1;
			ier=gauelm(m,1,wk,u,&det,iwk,m,&ifl);
			if(ier>0) {free(wk); free(v); free(iwk); return ier;}
		}
	}

	free(wk); free(v); free(iwk);
	return 142;
}
