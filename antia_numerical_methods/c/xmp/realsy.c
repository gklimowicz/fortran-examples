/*	Eigenvalues and eigenvectors of a real symmetric matrix */

#include <stdio.h>
#include <math.h>

int tred2(double *a, int n, int ia, double d[], double e[]);
int tql2(double *z, int n, int iz, double d[], double e[], double reps);
int trbak(double *a, int ia, int n, double *z, int iz, int nz);
double ran(double *seed);
int sturm(double e[], double d[], int n, int m1, int m2, double el[],
	double eu[], int *num, double reps);
int tinvit(double e[], double d[], int n, double el, double eu, double *ei,
	double ev[], double reps, int iflg, int *num);
int tridia(double e[], double d[], int n, int m1, int m2, double ei[],
	double eps1, double reps, double *ev, int iv);

main()
{
	int i,i1,j,n,m1,m2, id, iflg, ier,np,nmax;
	double hh, d[100],e[100], b[4][4],p,ei[20],eps,reps,p0;
	double a[4][4]={6,4,4,1, 4,6,1,4, 4,1,6,4, 1,4,4,6};

	id=4; n=4; np=0; reps=1.e-300;
	printf(" The matrix is :\n");
	for(i=0; i<n; ++i) {
		for(j=0; j<n; ++j) printf(" %e ",a[i][j]);
		printf("\n");
	}


/*	Reduce the matrix to tri-diagonal form */
	i=tred2(&a[0][0],n,id,d,e);
	printf(" ier = %d    n =  %d    diagonal elements :\n ",i,n);
	for(i=0; i<n; ++i) printf(" %e ",d[i]);
	printf(" \n off-diagonal elements : \n");
	for(i=0; i<n; ++i) printf(" %e ",e[i]);
	printf(" \n transformation matrix \n ");
	for(i=0; i<n; ++i) {
		for(j=0; j<n; ++j) printf(" %e ",a[i][j]);
		printf(" \n");
	}

		printf("type iflg =0/1  for QL algorithm/ Sturm sequence \n");
		scanf(" %d",  &iflg);

		if(iflg==0) {
/*	Use QL algorithm */
			reps=1.e-14;
			i=tql2(&a[0][0],n,id,d,e,reps);
			printf(" ier = %d \n   eigenvalue    eigenvector \n",i);
			for(i=0; i<n; ++i) {
				printf(" %e    ",d[i]);
				for(j=0; j<n; ++j) printf(" %e ",a[j][i]);
				printf(" \n");
			}
		}
		else {
/*	Find all eigenvalues and eigenvectors using Sturm sequence */
			m1=0; m2=n-1; eps=0.01;
			reps=1.e-14;
			i=tridia(e,d,n,m1,m2,ei,eps,reps,&b[0][0],id);
			np=m2-m1+1;
			printf(" ier = %d    m1 = %d    m2 = %d\n",i,m1,m2);
/*	Uncomment the following lines to print eigenvectors of tridiagonal matrix */

/*			printf("  eigenvalue    eigenvector \n");
			for(i=0; i<n; ++i) {
				printf(" %e    ",ei[i]);
				for(j=0; j<n; ++j) printf(" %e ",b[j][i]);
				printf(" \n");
			} */

			i=trbak(&a[0][0],id,n,&b[0][0],id,np);
			printf("  eigenvalue    eigenvector \n");
			for(i=0; i<n; ++i) {
				printf(" %e    ",ei[i]);
				for(j=0; j<n; ++j) printf(" %e ",b[j][i]);
				printf(" \n");
			}
		}
		

	return;
}




/*	To find specified eigenvalues of a real symmetric tridiagonal 
	matrix using bisection on Sturm sequence

	E : (input) Array of length N containing the off-diagonal elements
		of the tridiagonal matrix, E[i+1]=A(i,i+1)=A(i+1,i)
	D : (input) Array of length N containing the diagonal elements
		of the tridiagonal matrix, D[i]=A(i,i)
	N : (input) Order of the matrix
	M1 : (input) Serial number of lowest eigenvalue to be determined.
		The eigenvalues are sorted in increasing order
	M2 : (input) Serial number of highest eigenvalue to be determined.
		All eigenvalues from M1 to M2 are determined
	EL : (output) Array of length M2+1 containing the calculated
		lower limit on eigenvalues
	EU : (output) Array of length M2+1 containing the calculated
		upper limit on eigenvalues
		The ith eigenvalue is located in interval (EL[i],EU[i])
	NUM : (output) Number of times Sturm sequence was evaluated to locate
		the eigenvalues.
	REPS : (input) Relative accuracy to which eigenvalues are located by bisection
		
	Error status is returned by the value of the function STURM.
		0 value implies successful execution
		110, implies that M1<0 or M2>N-1, in which case no calculations are done

	Required functions : None
*/

#include <math.h>
#include <stdlib.h>

int sturm(double e[], double d[], int n, int m1, int m2, double el[],
	double eu[], int *num, double reps)

{
	int i,k,m,nl,nu,n1,n2,ifl,ns, nbis=70;
	double emin,emax,r,e1,q1,q2, xuf=-1.0e29, eps=1.e-30;
	double *wk;

	if(m1>m2) return 0;
	if(m1<0 || m2>=n) return 110;

	e[0]=0.0;
/*	Finding bounds on eigenvalues using Gerschgorin's theorem */
	emin=d[n-1]-fabs(e[n-1]);
	emax=d[n-1]+fabs(e[n-1]);
	wk=(double *) calloc((size_t) n, sizeof(double));
	wk[0]=0.0;
	for(i=n-2; i>=0; --i) {
		r=fabs(e[i])+fabs(e[i+1]);
		if(d[i]+r>emax) emax=d[i]+r;
		if(d[i]-r<emin) emin=d[i]-r;
		wk[i+1]=e[i+1]*e[i+1];
	}

/*	Initialise the limits to undefined values */
	for(i=m1; i<=m2; ++i) {el[i]=xuf; eu[i]=xuf;}
	nl=0; nu=n;
	*num=0;

/*	Loop for each eigenvalue */
	for(m=m1; m<=m2; ++m) {
/*	If the lower bound is undefined, then use EMIN */
		if(el[m]==xuf) {el[m]=emin; n1=nl;}

		if(eu[m]==xuf) {
/*	If upper bound is undefined, use the bound for some higher eigenvalue
	and if none exists, then use EMAX */
			for(i=m+1; i<=m2; ++i) {
				if(eu[i] != xuf) {
					eu[m]=eu[i];
					n2=i;
					break;
				}
			}
			if(eu[m]==xuf) {eu[m]=emax; n2=nu;}
		}

		ifl=0;
/*	Loop for bisection */
		for(i=1; i<=nbis; ++i) {
			e1=(el[m]+eu[m])/2.;
			if(e1==el[m] || e1==eu[m]) break;
			*num=(*num)+1;

/*	Count the number of sign changes in the Sturm sequence */
			ns=0;
			q1=d[0]-e1;
			if(q1<0.0) ns=1;
			for(k=1; k<n; ++k) {
				q2=d[k]-e1;
				if(q1 != 0.0) q2=q2-wk[k]/q1;
				else q2=q2-fabs(e[k])/eps;
				if(q2<0.0) ns=ns+1;
				q1=q2;
			}

/*	If the bounds are two consecutive real number on the machine, then quit */
			if(e1 == el[m] || e1==eu[m]) break;
/*	Update the bounds */
			if(ns>=m1+1 && ns<=m2+1) {
				if(e1<eu[ns-1]) eu[ns-1]=e1;
			}
			if(ns>=m1 && ns<m2+1) {
				if(e1>el[ns]) el[ns]=e1;
			}
			if(ns<m+1) {
				if(ns>n1) n1=ns;
				if(e1>el[m]) el[m]=e1;
			}
			else {
				if(ns<n2) n2=ns;
				if(e1<eu[m]) eu[m]=e1;
			}

			if(n1==m && n2==m+1) {
/*	The eigenvalue is isolated */
				if(ifl>3 && fabs(eu[m]-el[m])) goto next;
				if(m==m1) ifl=ifl+1;
				else if(el[m] != eu[m-1]) ifl=ifl+1;
			}
		}

/*	If the eigenvalue cannot be isolated, then set the same bounds
	for all of them */
		for(k=m+1; k<n2; ++k) {el[k]=el[m]; eu[k]=eu[m];}

next:	n1=m+1; if(m+2>n2) n2=m+2;
	}
	free(wk);
	return 0;
}



/*	To find a specified eigenvalue and eigenvector of a real symmetric
	tridiagonal matrix using inverse iteration

	E : (input) Array of length N containing the off-diagonal elements
		of the tridiagonal matrix, E[i+1]=A(i,i+1)=A(i+1,i)
	D : (input) Array of length N containing the diagonal elements
		of the tridiagonal matrix, D[i]=A(i,i)
	N : (input) Order of the matrix
	EL : (input) lower limit on eigenvalue
	EU : (input) upper limit on eigenvalue
	EI : (output) The calculated eigenvalue
	EV : (input/output) Array of length N containing the calculated
		eigenvector. If IFLG!=0 then EV should contain the previous
		eigenvector, so that initial guess to EV is orthogonal to it
	REPS : (input) Relative accuracy to which eigenvalues or eigenvectors
		are calculated
	IFLG : (input) Integer variable used as a flag to denote how initial
		approximation to eigenvector is obtained. If IFLG=0 then
		initial guess to eigenvector is generated using random numbers.
		Otherwise, initial guess is chosen to be orthogonal to EV
	NUM : (output) Number of iterations required by the function
		
	Error status is returned by the value of the function TINVIT.
		0 value implies successful execution
		144 implies that inverse iteration failed to converge

	Function RAN1(SEED) is required to generate random numbers.
	If a different function is used, the seed should be changed
	appropriately.

	Required functions : RAN1
*/

#include <math.h>
#include <stdlib.h>

double ran1(double *seed);

int tinvit(double e[], double d[], int n, double el, double eu, double *ei,
	double ev[], double reps, int iflg, int *num)

{
	int i,j,k,km, nit=100;
	double s,evs,el1,eu1,t,r,r1,s1,s2,dei,ukm,dvs,rep,hi;
	double *u, *b;
	static double seed=12345;   /*	The seed for random number generator */

	e[0]=0.0;
	u=(double *) calloc((size_t) n,sizeof(double));
	b=(double *) calloc((size_t) (4*n),sizeof(double));

	if(iflg==0) {
/*	Choose a random initial guess for eigenvector */
		for(i=0; i<n; ++i) {
			u[i]=ran1(&seed);
			ev[i]=u[i];
		}
	}

	else {
/*	Choose the initial vector to be orthogonal to EV */
		s=0.0; evs=0.0;
		for(i=0; i<n; ++i) {
			u[i]=ran1(&seed);
			evs=evs+ev[i]*ev[i];
			s=s+u[i]*ev[i];
		}
		s=s/evs;
		for(i=0; i<n; ++i) {
			ev[i]=u[i]-ev[i]*s;
			u[i]=ev[i];
		}
	}

	*ei=(el+eu)/2.;
/*	Stretch the bounds to allow for roundoff errors */
	el1=el-reps*fabs(el);
	eu1=eu+reps*fabs(eu);

/*	Loop for inverse iteration */
	for(j=1; j<=nit; ++j) {
		*num=j;

/*	Set up the matrix A-pI */
		for(i=0; i<n; ++i) {
			b[i*4]=e[i];
			b[1+i*4]=d[i]-(*ei);
			if(i<n) b[2+i*4]=e[i+1];
			b[3+i*4]=0.0;
		}
		b[4*n-2]=0.0;

/*	Gaussian elimination with partial pivoting */
		for(k=0; k<n-1; ++k) {
			if(fabs(b[1+k*4])<fabs(b[4*k+4])) {
				for(i=0; i<3; ++i) {
					t=b[i+4*k+4]; b[i+4*k+4]=b[i+1+4*k]; b[i+1+4*k]=t;
				}
				t=ev[k+1]; ev[k+1]=ev[k]; ev[k]=t;
			}

/*	If a pivot is zero make it nonzero */
			if(b[1+4*k]==0.0) b[1+4*k]=reps;
			r=-b[4*k+4]/b[1+4*k];
			b[4*k+5]=b[4*k+5]+r*b[2+4*k];
			b[4*k+6]=b[4*k+6]+r*b[3+4*k];
			ev[k+1]=ev[k+1]+r*ev[k];
		}
		if(b[4*n-3]==0.0) b[4*n-3]=reps;

/*	Back-substitution */
		ev[n-1]=ev[n-1]/b[4*n-3];
		ev[n-2]=(ev[n-2]-ev[n-1]*b[4*n-6])/b[4*n-7];
		for(k=n-3; k>=0; --k) ev[k]=(ev[k]-ev[k+1]*b[2+4*k]-ev[k+2]*b[3+4*k])/b[1+4*k];

/*	Calculating the Rayleigh quotient (DEI) */
		r1=0.0; km=0; s1=0.0; s2=0.0;
		for(k=0; k<n; ++k) {
			s1=s1+ev[k]*ev[k];
			s2=s2+ev[k]*u[k];
			if(r1<fabs(ev[k])) {r1=fabs(ev[k]); km=k;}
		}

		dei=s2/s1;
		ukm=ev[km];
		dvs=0.0;
/*	Normalising the eigenvector */
		for(k=0; k<n; ++k) {
			ev[k]=ev[k]/ukm;
			dvs=dvs+fabs(ev[k]-u[k]);
		}

/*	The convergence check */
		rep=reps*fabs(*ei);
		if(fabs(dei)<rep || r1*rep>1.0 || dvs<reps) {free(b); free(u); return 0;}
		for(i=0; i<n; ++i) u[i]=ev[i];
		if(*ei+dei>el1 && (*ei)+dei<eu1) *ei=(*ei)+dei;
		else {
/*	If the new shift is outside the specified bounds, then modify it */
			hi=fabs(*ei-el); if(hi>fabs(*ei-eu)) hi=fabs(*ei-eu);
			s1=hi*0.5; if(dei<0.0) s1=-s1;
			*ei=(*ei)+s1;
		}
	}

	free(b); free(u);
	return 144;
}



/*	To find eigenvalues and eigenvectors of Z T Z^T using QL algorithm
	where T is a symmetric tridiagonal matrix and Z is an orthogonal matrix.
	If Z is the transformation matrix to reduce original matrix to
	tridiagonal matrix, it will calculate the eigenvectors of original matrix

	Z : (input/output) Array of length IZ*N which should contain
		the transformation matrix required to reduce original real
		symmetric matrix to tridiagonal form. To find eigenvectors
		of a symmetric tridiagonal matrix, set Z to unit matrix.
		After execution Z will contain the eigenvector of the original
		matrix Z T Z^T. Z[i][j] should contain the ith component of
		jth eigenvector
	N : (input) Order of matrix
	IZ : (input) The second dimension of array Z as declared in the
		calling function. (IZ>=N)
	D : (input/output) Array of length N, containing the diagonal
		elements of the tridiagonal matrix, D[i]=T[i][i]
		After execution it will contain the eigenvalues of the matrix
	E : (input/output) Array of length N containing the off-diagonal
		elements of the tridiagonal matrix, E[i]=T[i][i+1]=T[i+1][i]
		It is used as scratch space and its contents will be destroyed
		during execution.
	REPS : (input) Required tolerance, it should be of order of machine
		accuracy
		
	Error status is returned by the value of the function TQL2.
		0 value implies successful execution
		108 implies that N<=1 or N>IZ, in which case no
			calculations are performed
		143 implies that the QL algorithm failed to converge
			for some eigenvalue, the calculations are abandoned

	Required functions : None
*/

#include <math.h>

int tql2(double *z, int n, int iz, double d[], double e[], double reps)

{
	int i,j,k,l,m,it, nit=30;
	double b,f,g,h,p,r,c,s;

	if(n<=0 || n>iz) return 108;

	for(i=1; i<n; ++i) e[i-1]=e[i];
	e[n-1]=0.0;
	b=0.0; f=0.0;

	for(l=0; l<n; ++l) {
		h=reps*(fabs(d[l])+fabs(e[l]));
		if(b<h) b=h;
/*	Look for small off-diagonal elements */
		for(i=l; i<n; ++i) {
			m=i;
			if(fabs(e[m])<=b) break;
		}
		if(m==l) {
/*	one eigenvalue is isolated */
			d[l]=d[l]+f;
			goto nextit;
		}

/*	Loop for QL transformation  */
		for(it=1; it<=nit; ++it) {
/*	Find shift */
			g=d[l];
			p=(d[l+1]-g)/(2.0*e[l]);
			r=sqrt(p*p+1.0);
			if(p<0.0) r=-r;
			d[l]=e[l]/(p+r);
			h=g-d[l];
			for(i=l+1; i<n; ++i) d[i]=d[i]-h;
			f=f+h;

/*	The QL transformation */
			p=d[m];
			c=1.0; s=0.0;
/*	Given's rotations */
			for(i=m-1; i>=l; --i) {
				g=c*e[i];
				h=c*p;
				if(fabs(p)>=fabs(e[i])) {
					c=e[i]/p;
					r=sqrt(c*c+1.0);
					e[i+1]=s*p*r;
					s=c/r;
					c=1.0/r;
				}
				else {
					c=p/e[i];
					r=sqrt(c*c+1.0);
					e[i+1]=s*e[i]*r;
					s=1.0/r;
					c=c/r;
				}
				p=c*d[i]-s*g;
				d[i+1]=h+s*(c*g+s*d[i]);

/*	Transforming the eigenvectors */
				for(k=0; k<n; ++k) {
					h=z[i+1+k*iz];
					z[i+1+k*iz]=s*z[i+k*iz]+c*h;
					z[i+k*iz]=c*z[i+k*iz]-s*h;
				}
			}
			e[l]=s*p;
			d[l]=c*p;
			if(fabs(e[l])<=b) {
/*	One eigenvalue is isolated */
				d[l]=d[l]+f;
				goto nextit;
			}
		}
/*	QL iteration fails to converge */
		return 143;

nextit: b=b;
	}

/*	Sort eigenvalues in ascending order by straight selection */
	for(i=0; i<n-1; ++i) {
		k=i;
		p=d[i];
		for(j=i+1; j<n; ++j) {
			if(d[j]<p) {k=j; p=d[j];}
		}
		if(k != i) {
/*	exchange the eigenvalues and eigenvectors */
			d[k]=d[i];
			d[i]=p;
			for(j=0; j<n; ++j) {
				p=z[i+j*iz]; z[i+j*iz]=z[k+j*iz]; z[k+j*iz]=p;
			}
		}
	}
	return 0;
}



/*	To perform back-transformation on eigenvectors of reduced tridiagonal matrix
	to obtain eigenvectors of original real symmetric matrix reduced by TRED2
	This function is not required if eigenvectors are computed using TQL2

	A : (input/output) Array of length IA*N containing the transformation
		matrix. The last column of this array is used as scratch
		space and will be destroyed during execution
	IA : (input) The second dimension of array A as declared in the
		calling function (IA>=N)
	N : (input) Order of the matrix
	Z : (input/output) Array of length IZ*N containing the
		eigenvectors of tridiagonal matrix. After execution it will
		be overwritten by eigenvectors of original matrix
	IZ : (input) Second dimension of array Z as declared in the calling
		function (IZ>=NZ)
	NZ : (input) Number of eigenvectors to be transformed

	The returned value is always zero.

	Required functions : None
*/

#include <math.h>

int trbak(double *a, int ia, int n, double *z, int iz, int nz)

{
	int i,j,k;
	double s;

/*	Loop on eigenvectors */
	for(i=0; i<nz; ++i) {
/*	The matrix multiplication */
		for(j=0; j<n; ++j) {
			s=0.0;
			for(k=0; k<n-1; ++k) s=s+a[k+j*ia]*z[i+k*iz];
/*	To take care of the last column of A, which is overwritten */
			if(j==n-1) s=s+z[i+(n-1)*iz];
			a[n-1+j*ia]=s;
		}

		for(j=0; j<n; ++j) z[i+j*iz]=a[n-1+j*ia];
	}
	return 0;
}



/*	Reduction of real symmetric matrix to tridiagonal form using
	Householder's method

	A : (input/output) Array of length IA*N containing the matrix
		elements. After execution it will be overwritten by the
		transformation matrix
	N : (input) Order of matrix
	IA : (input) The second dimension of A as specified in the calling function
	D : (output) Diagonal elements of the transformed tridiagonal matrix
		D[I] would contain A[I][I].
	E : (output) Off-diagonal elements of the transformed tridiagonal matrix
		E[I+1] would contain A[I][I+1]=A[I+1][I]
		
	Error status is returned by the value of the function TRED2.
		0 value implies successful execution
		107 implies that N<=1 or N>IA, in which case no calculations are done

	Required functions : None
*/

#include <math.h>

int tred2(double *a, int n, int ia, double d[], double e[])

{
	int i,j,k;
	double f,g,h,hh, reps=1.e-300;

	if(n<=1 || n>ia) return 107;

	for(i=n-1; i>=1; --i) {
		f=a[i-1+i*ia];
		g=0.0;
		for(k=0; k<=i-2; ++k) g=g+a[k+i*ia]*a[k+i*ia];
		h=g+f*f;
		if(g<= reps) {
/*	Skip the transformation */
			e[i]=f;
			h=0.0;
		}

		else {
			g=sqrt(h);
			if(f>0.0) g=-g;
			e[i]=g;
			h=h-f*g;
			a[i-1+i*ia]=f-g;
			f=0.0;

			for(j=0; j<=i-1; ++j) {
/*	Elements of u_i/H_i */
				a[i+j*ia]=a[j+i*ia]/h;
				g=0.0;
/*	Form elements of A_iu_i */
				for(k=0; k<=j; ++k) g=g+a[k+j*ia]*a[k+i*ia];
				for(k=j+1; k<=i-1; ++k) g=g+a[j+k*ia]*a[k+i*ia];
/*	Components of p_i */
				e[j]=g/h;
				f=f+g*a[i+j*ia];
			}

/*	calculate u_i^Tp_i/2H_i */
			hh=0.5*f/h;
			for(j=0; j<=i-1; ++j) {
				f=a[j+i*ia];
				g=e[j]-hh*f;
/*	Elements of q_i */
				e[j]=g;
				for(k=0; k<=j; ++k) a[k+j*ia]=a[k+j*ia]-f*e[k]-g*a[k+i*ia];
			}
		}
		d[i]=h;
	}

	d[0]=0.0; e[0]=0.0;
/*	accumulation of transformation matrix Q */
	for(i=0; i<n; ++i) {
		if(d[i] != 0.0) {
			for(j=0; j<=i-1; ++j) {
				g=0.0;
				for(k=0; k<=i-1; ++k) g=g+a[k+i*ia]*a[j+k*ia];
				for(k=0; k<=i-1; ++k) a[j+k*ia]=a[j+k*ia]-g*a[i+k*ia];
			}
		}
		d[i]=a[i+i*ia];
		a[i+i*ia]=1.0;
		for(j=0; j<=i-1; ++j) {a[j+i*ia]=0.0; a[i+j*ia]=0.0;}
	}
	return 0;
}



/*	To find specified eigenvalues and eigenvectors of a real symmetric
	tridiagonal matrix using Sturm sequence coupled with inverse iteration

	E : (input) Array of length N containing the off-diagonal elements
		of the tridiagonal matrix, E[i+1]=A(i,i+1)=A(i+1,i)
	D : (input) Array of length N containing the diagonal elements
		of the tridiagonal matrix, D[i]=A(i,i)
	N : (input) Order of the matrix
	M1 : (input) Serial number of lowest eigenvalue to be determined.
		The eigenvalues are sorted in increasing order
	M2 : (input) Serial number of highest eigenvalue to be determined.
		All eigenvalues from M1 to M2 are determined
	EI : (output) Array of length M2-M1+1 containing the calculated
		eigenvalues
	EPS1 : (input) Relative accuracy to which eigenvalues are located by bisection
		before using inverse iteration. If the inverse iteration
		does not converge to nearest eigenvalues, EPS1 can be reduced.
		A value of 0.1-0.01 times typical eigenvalue is generally sufficient.
	REPS : (input) Desired relative accuracy in eigenvalues and eigenvectors
	EV : (output) Array of length IV*N containing the
		eigenvectors. EV[i][j] should contain the ith component of
		the jth eigenvector
	IV : (input) The second dimension of array EV as declared in the
		calling function (IV>=M2-M1+1)
		
	Error status is returned by the value of the function TRIDIA.
		0 value implies successful execution
		109 implies that N<=1 or M2-M1+1>IV or M1<0 or M2>=N
		Other values may be set by TINVIT, only the last value is
		returned.

	Required functions : STURM, TINVIT, RAN1
*/

#include <math.h>
#include <stdlib.h>

int sturm(double e[], double d[], int n, int m1, int m2, double el[],
	double eu[], int *num, double reps);
int tinvit(double e[], double d[], int n, double el, double eu, double *ei,
	double ev[], double reps, int iflg, int *num);

int tridia(double e[], double d[], int n, int m1, int m2, double ei[],
	double eps1, double reps, double *ev, int iv)

{
	int i,j,i1,iflg,ier,num;
	double *wl,*wu,*v;

	if(n<=1 || m2-m1+1>iv || m1<0 || m2>=n) return 109;

	if(m1>m2) return 0;

	wl=(double *) calloc((size_t) n,sizeof(double));
	wu=(double *) calloc((size_t) n,sizeof(double));
	v=(double *) calloc((size_t) n,sizeof(double));
/*	Locate the eigenvalues */
	ier=sturm(e,d,n,m1,m2,wl,wu,&num,eps1);
	ier=0;

/*	Loop for finding individual eigenvalues and eigenvectors */
	for(i=m1; i<=m2; ++i) {
		iflg=0;
		i1=i-m1;

		if(i>m1) {
			if(fabs(wl[i]-wl[i-1])<3.*reps*fabs(wl[i])) {
/*	Set the flag for close eigenvalues */
				iflg=1;
				for(j=0; j<n; ++j) v[j]=ev[i1-1+j*iv];
			}
		}
		j=tinvit(e,d,n,wl[i],wu[i],&ei[i1],v,reps,iflg,&num);
		if(j>0) ier=j;
		for(j=0; j<n; ++j) ev[i1+j*iv]=v[j];
	}
	free(v); free(wu); free(wl);
	return ier;
}



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
	
