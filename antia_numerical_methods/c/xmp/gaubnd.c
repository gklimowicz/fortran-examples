/*	Solution of linear algebraic equations with a band matrix using
	Gaussian elimination */

#include <stdio.h>
#include <math.h>

int gaubnd(int n, int kb, int num, double *a, double *x, double *det,
		int *idet, int inc[], int lj, int *iflg);
double ranf(int *iseed);

main()
{
	int i,i1,j,k,n, num, iflg, is, lj, inc[20], idet;
	double r1, det, x[2][20], wk[60], am[20][20], bm[20][20];

	num=1;
	lj=20;
	is=-2;
	for(i1=0; i1<99; ++i1) {
		printf("type n = no. of eqs,  k = band width     (quits when n<=0)\n");
		scanf("%d %d", &n, &k);
		if(n<=0) return 0;

/*	Generate a matrix using random numbers and RHS assuming solution =1 */
		for(i=0; i<n; ++i) {
			r1=0.0;
			for(j=0; j<2*k+1; ++j) {
				am[j][i]=0.0;
				if(j>=k-i && j<k+n-i) {
					am[j][i]=ranf(&is);
					r1=r1+am[j][i];
				}
			}
			x[0][i]=r1;
        }
		iflg=0;
		i=gaubnd(n,k,num,&am[0][0],&x[0][0],&det,&idet,inc,lj,&iflg);
		printf(" ier = %d  n = %d  k = %d   det = %e %d \n", i,n,k,det,idet);
		printf(" solution : ");
		for(j=0; j<n; j++) printf(" %f ",x[0][j]);
		printf("\n");
	}
	return;
}

 

/*     Solution of a system of linear equations using Gaussian elimination
     	for a band matrix

	N : (input) Number of equations to be solved
	KB : (input) Bandwidth of matrix A[i][j]=0 if abs(I-J)>KB
	NUM : (input) Number of different sets (each with N equations) of
		equations to be solved
	A : (input/output) The matrix of coefficients of size LJ*(3*KB+1)
		A[J-I+KB][I] is the coefficient of x_J in Ith equation
		at output it will contain the triangular decomposition
	X : (input/output) The matrix containing right hand sides (size LJ*NUM)
		X[j][i] is the ith element of jth right hand side
		at output it will contain the solutions
	DET, IDET : (output) The determinant of the matrix = DET*2**IDET
	INC : (output) Integer array of length N containing information about
		interchanges performed during elimination
	LJ : (input) Second dimension of arrays A and X in calling function
	IFLG : (input) Integer variable used as a flag to specify the type
		of computation required
     	If IFLG=-1, both elimination and solution are calculated
     		without pivoting and IFLG is set to 2
	If IFLG=0, both elimination and solution are computed
     		with partial pivoting and IFLG is set to 2
     	If IFLG=1, only elimination is done with pivoting and IFLG is set to 2
     	If IFLG>=2 only solution is calculated, the triangular
     		decomposition should have been calculated earlier
		
	Error status is returned by the value of the function GAUBND.
		0 value implies successful execution
		104 implies N<=0 or N>LJ or KB>N
		124 implies some pivot turned out to be zero and hence
	     		matrix must be nearly singular

	Required functions : None
*/

#include <math.h>
#include <stdlib.h>

int gaubnd(int n, int kb, int num, double *a, double *x, double *det,
		int *idet, int inc[], int lj, int *iflg)

{
	int i,j,k,km,l,kb1,m1,m2;
	double r1, t1;
	double *wk;

	if(n<=0 || n>lj || kb>n) return 104;

	kb1=kb+1;
	if(*iflg < 2) {
/*     Perform elimination */
		for(i=0; i<n ; ++i) {
			for(j=2*kb+1; j<=3*kb; ++j) a[j*lj+i] = 0.0;
		}

		*det=1.0; *idet=0;
		wk=(double *) calloc((size_t) (3*kb+1), sizeof(double));
		for(k=0; k<n-1; ++k) {
/*     Find the maximum element in the Kth column  */
			r1=0.0; km=k;
			if(*iflg >= 0) {
				m1=k+kb; if(n-1 < m1) m1=n-1;
				for(l=k; l<=m1; ++l) {
					if(fabs(a[l+(k-l+kb)*lj]) > r1) {
						r1=fabs(a[l+(k-l+kb)*lj]);
						km=l;
					}
				}
			}

			inc[k]=km;
			if(km != k) {
/*     Interchange the rows if needed  */
				m1=2*kb+k; if(n-1 < m1) m1=n-1;
				for(l=k; l<=m1; ++l) wk[l-k]=a[k+(l-k+kb)*lj];
				for(l=k; l<=m1; ++l) {
					a[k+(l-k+kb)*lj]=a[km+(l-km+kb)*lj];
					a[km+(l-km+kb)*lj]=wk[l-k];
				}
				*det= -(*det);
			}

			*det = (*det)*a[kb*lj+k];
			if( a[kb*lj+k] == 0.0) {free(wk); return 124;}
			if(*det != 0.0) {
/*     Scale the value of the determinant   */
				while(fabs(*det) > 32.0) {
					*det = (*det)*0.03125e0; *idet = *idet + 5;
				}

				while(fabs(*det) < 0.03125e0) {
					*det = (*det)*32.0; *idet = *idet - 5;
				}
			}

			m1=k+kb; if(n-1 < m1) m1=n-1;
			for(l=k+1; l<= m1; ++l) {
				a[l+(k-l+kb)*lj] = a[l+(k-l+kb)*lj]/a[kb*lj+k];
				m2=k+2*kb; if(n-1 < m2) m2=n-1;
				for(i=k+1; i<=m2; ++i)
					a[l+(i-l+kb)*lj]=a[l+(i-l+kb)*lj]-a[l+(k-l+kb)*lj]*a[k+(i-k+kb)*lj];
			}
		}

		free(wk);
		*det = (*det)*a[(n-1)+kb*lj];
		inc[n-1]=n-1;
		if(a[(n-1)+kb*lj] == 0.0) return 124;

		if(*iflg==1) {*iflg=2; return 0;}
		*iflg=2;
	}
		
/*     Solution for the NUM different right-hand sides */
	for(j=0; j<num; ++j) {
/*     Forward substitution  */
		for(k=0; k<n-1; ++k) {
			if(k != inc[k]) {
				t1= x[j*lj+k];
				x[j*lj+k] = x[j*lj+inc[k]];
				x[j*lj+inc[k]] = t1;
			}
			m1=k+kb; if(n-1 < m1) m1=n-1;
			for(l=k+1; l<=m1; ++l)
				x[j*lj+l]=x[j*lj+l]-a[l+(k-l+kb)*lj]*x[j*lj+k];
		}

/*     back-substitution  */
		x[j*lj+n-1] = x[j*lj+n-1]/a[(n-1)+kb*lj];
		for(k=n-2; k>=0; --k) {
			m1=k+2*kb; if(n-1 < m1) m1=n-1;
			for(l=m1; l>=k+1; --l)
				x[j*lj+k]=x[j*lj+k]-x[j*lj+l]*a[k+(l-k+kb)*lj];
			x[j*lj+k] = x[j*lj+k]/a[kb*lj+k];
		}
	}
	return 0;
}





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
