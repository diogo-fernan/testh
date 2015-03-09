
#include "func.h"
#include "io.h"

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>

#define EPS1 	(double) 1.0e-6
#define JMAX 	40

double func_qsimp (
	dist 	d, 
	double 	a, 
	double 	b)
// Returns the integral of the function func from a to b.
// The parameters EPS can be set to the desired fractional accuracy and 
// JMAX so that 2 to the power JMAX-1 is the maximum allowed number of steps
// Integration is performed by Simpson’s rule.
{
	double trapzd (dist d, double a, double b, int n); 

	int j;
	double s, st, ost=0.0, os=0.0;

	for (j=1; j<=JMAX; j++) { 
		st = trapzd (d, a, b, j);
		s = (4.0*st-ost) / 3.0;
		if (j > 5)
			if (fabs (s-os) < EPS1*fabs (os) || (s == 0.0 && os == 0.0)) 
				return s;
		os = s;
		ost = st;
	}
	io_PrintError (ERR, "no convergence in"
			" func_qsimp"); 
	return 0.0;
}

double trapzd (
	dist 	d, 
	double 	a, 
	double 	b, 
	int 	n)
// This routine computes the nth stage of refinement of an 
// extended trapezoidal rule. func is input￼as a pointer to the function to be 
// integrated between limits a and b, also input. When called with
// n=1, the routine returns the crudest estimate of \int_a^b f (x)dx. 
// Subsequent calls with n=2,3,... (in that sequential order) will 
// improve the accuracy by adding 2^(n-2) additional interior points. 
{
	double x, tnm, sum, del; 
	static double s;
	int it, j;
	if (n == 1) {
		s = 0.5*(b-a) * (dist_PDFProb (d, a) + dist_PDFProb (d, b));
		return s;
	} else {
		for (it=1, j=1; j<n-1; j++) 
			it <<= 1;
		tnm = it;
		del = (b-a) / tnm; 
		x = a + 0.5*del;
		for (sum=0.0, j=1; j<=it; j++, x+=del) 
			sum += dist_PDFProb (d, x); 
		s = 0.5 * (s + (b-a) * sum/tnm);
		return s;
	} 
}

#define ITMAX 	100 
#define EPS2 	3.0e-7 
#define FPMIN 	1.0e-30

double func_gammq (
	double a, 
	double x)
// Returns the incomplete gamma function Q(a, x) = 1 − P (a, x)
{
	if (x < 0.0 || a <= 0.0) 
		io_PrintErr (ERR, "invalid parameters in"
			" func_gammq"); 

	void gcf (double *gammcf, double a, double x, double *gln); 
	void gser (double *gamser, double a, double x, double *gln);

	double gamser, gammcf, gln;

	if (x < a + 1.0) { 
		gser (&gamser, a, x, &gln); 
		return 1.0 - gamser;
	} else { 
		gcf (&gammcf, a, x, &gln); 
		return gammcf;
	} 
}

float gammln (
	double xx) 
// Returns the value ln[Γ(xx)] for xx > 0
{
	int j;
	double x, y, tmp, ser;
	
	static double cof[6]={76.18009172947146, -86.50532032941677,
	24.01409824083091, -1.231739572450155,
	0.1208650973866179e-2, -0.5395239384953e-5}; 
	
	y = x = xx;
	tmp = x + 5.5;
	tmp -= (x+0.5) * log (tmp); 
	ser = 1.000000000190015;
	for (j=0; j<=5; j++) 
		ser += cof[j] / ++y; 
	return -tmp + log (2.5066282746310005*ser / x);
}

void gcf (
	double *gammcf, 
	double a, 
	double x, 
	double *gln)
// Returns the incomplete gamma function Q(a, x) evaluated by 
// its continued fraction representation as gammcf 
// Also returns lnΓ(a) as gln
{
	int i;
	double an, b, c, d, del, h;
	
	*gln = gammln (a);
	b = x + 1.0 - a;
	c = 1.0 / FPMIN;
	d = 1.0 / b;
	h = d;

	for (i=1; i<=ITMAX; i++) {
		an = -i * (i-a);
		b += 2.0;
		d = an * d + b;
		if (fabs (d) < FPMIN)
			d = FPMIN;
		c = b + an / c;
		if (fabs (c) < FPMIN)
			c = FPMIN;
		d = 1.0 / d;
		del = d * c;
		h *= del;
		if (fabs (del-1.0) < EPS2) 
			break;
	}

	if (i > ITMAX) 
		io_PrintError (ERR, "a too large, ITMAX too small in"
			" gcf"); 
	*gammcf = exp (-x + a * log (x) - *gln) * h; 
}

void gser (
	double *gamser, 
	double a, 
	double x, 
	double *gln)
// Returns the incomplete gamma function P (a, x) evaluated by 
// its series representation as gamser
// Also returns lnΓ(a) as gln
{
	if (x <= 0.0)
		io_PrintErr (ERR, "invalid parameters in"
			" gser");  

	int n;
	double sum, del, ap;

	*gln = gammln (a); 

	ap = a;
	del = sum = 1.0 / a;
	for (n=1; n<=ITMAX; n++) {
		++ap;
		del *= x / ap;
		sum += del;
		if (fabs (del) < fabs (sum) * EPS2) {
			*gamser = sum * exp (-x + a * log (x) - *gln);
			return; 
		}
	}
	io_PrintError (ERR, "a too large, ITMAX too small in"
		" gser"); 
}

long func_GCDR (long a, long b)
{
	if (a == 0) return b;
	return func_GCDR (b % a, a);
}
