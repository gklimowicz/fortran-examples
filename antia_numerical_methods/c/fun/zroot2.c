/*	Complex zeros of a given function using Muller iteration with deflation
	For use with MULER2 which requires function in form CF*2**IDET

	N : (input) Number of zeros to be determined
	CX : (input) Array of length N*2 containing the initial guesses
		for the complex zeros
	CZERO : (input/output) Array of length (N+NZ)*2 containing the
		computed values of the complex zeros. The first NZ zeros which are
		already known must be supplied as input while other zeros will be added
	NZ : (input/output) Number of known zeros. At input it should
		contain the number already known which are included in
		array CZERO. This number will be incremented as more zeros
		are determined successfully.
	REPS : (input) Required relative accuracy
	AEPS : (input) Required absolute accuracy
		The estimated error should be less than max(AEPS, REPS*fabs(X))
	RMAX : (input) Maximum magnitude of zeros. Iteration will be terminated
		when ABS(CZ) > RMAX
	CF : (input) Name of the function routine to calculate the function
		Function CF(Z,F,JF) must be supplied by the user, where F and Z
		are arrays of length 2 containing the complex values. F*2**JF should
		be the calculated value of the complex function at Z.
		
	Error status is returned by the value of the function ZROOT2.
		0 value implies successful execution
		41 implies that iteration did not converge to satisfactory accuracy for
			at least one zero, but value is acceptable at lower accuracy.
		429 implies that iteration failed to find at least one zero.
			The number of zeros found can be checked from the value of NZ.

	Required functions : MULER2, CF, CABS, CDIV, CSQRT
*/

#include <math.h>

void muler2(double cx1[], double cx2[], double cx3[], double reps, double aeps,
	int *ier, double cf[], double cx[], int ix, int nz, double czero[][2], double rmax);

int zroot2(int n, double cx[][2], double czero[][2], int *nz, double reps,
	double aeps, double rmax, void cf(double * , double * , int * ))

{
	int i,ier,ier1,ix;
	double cdx[2],cx1[2],cx2[2],cx3[2],cf0[2],cx0[2];

	if(n<=0) return 0;
	if(*nz<0) *nz=0;
	ier=0;

	for(i=0; i<n; ++i) {
/*	The starting values for Muller's iteration */
		cx3[0]=cx[i][0]; cx3[1]=cx[i][1];
		cdx[0]=0.01*cx3[0]; cdx[1]=0.01*cx3[1];
		if(i<n-1 && cdx[0]==0.0 && cdx[1]==0.0) {
			cdx[0]=0.1*(cx[i+1][0]-cx[i][0]);
			cdx[1]=0.1*(cx[i+1][1]-cx[i][1]);
		}
		if(i>0 && cdx[0]==0.0 && cdx[1]==0.0) {
			cdx[0]=0.1*(cx[i-1][0]-cx[i][0]);
			cdx[1]=0.1*(cx[i-1][1]-cx[i][1]);
		}
		if(cdx[0]==0.0 && cdx[1]==0.0) {cdx[0]=1.0e-3; cdx[1]=0.0;}
/*	These values may need to be changed in some cases */
		cx2[0]=cx3[0]+cdx[0]; cx2[1]=cx3[1]+cdx[1];
		cx1[0]=cx3[0]-cdx[0]; cx1[1]=cx3[1]-cdx[1];

/*	Find the next zero */
		ier1=0;
		do {
			muler2(cx1,cx2,cx3,reps,aeps,&ier1,cf0,cx0,ix,*nz,czero,rmax);
			if(ier1>=0) break;
			cf(cx0,cf0,&ix);
		} while(ier1<0);

		if(ier1<100) {
/*	The zero is accepted */
			czero[*nz][0]=cx3[0]; czero[*nz][1]=cx3[1];
			*nz=(*nz)+1;
			if(ier==0 && ier1 != 0) ier=41;
		}
		else ier=429;
	}
	return ier;
}
