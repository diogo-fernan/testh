
#include "stat.h"
#include "dist.h"
#include "rng.h"
#include "stat.h"
#include "util.h"
#include "io.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>


int dist_CheckDist (
	dist d)
{
	int i, f;
	for (i=TestH_Gauss, f=0; i<=TestH_Exp; i++) {
		if (d == i) {
			f = 1;
			break;
		}
	}
	if (!f)
		return ERR;
	return OK;
}

double dist_PDFProb (
	dist 	d,
	double 	x)
{
	switch (d) {
		case TestH_Unif:
			return dist_UnifPDF (0.0, 1.0);
			break;
		case TestH_Gauss:
			return dist_GaussPDF (x, 0.0, 1.0);
			break;
		case TestH_Pareto:
			return dist_ParetoPDF (x, 1.0, 1.0);
			break;
		case TestH_Exp:
			return dist_ExponentialPDF (x, 1.0);
			break;
		default:
			return 0.0;
			break;
	}
}
double dist_CDFProb (
	dist 	d,
	double 	x)
{
	switch (d) {
		case TestH_Unif:
			return dist_UnifCDF (x, 0.0, 1.0);
			break;
		case TestH_Gauss:
			return dist_GaussCDF (x, 0.0, 1.0);
			break;
		case TestH_Pareto:
			return dist_ParetoCDF (x, 1.0, 1.0);
			break;
		case TestH_Exp:
			return dist_ExponentialCDF (x, 1.0);
			break;
		default:
			return 0.0;
			break;
	}
}

double dist_UnifPDF (
	double a,
	double b)
{
	return 1.0 / (b - a);
}
double dist_GaussPDF (
	double x,
	double mu,
	double sigma)
{
	double y = 1 / (sigma*sqrt (2*PI)) * 
		pow (NEPER, -((x-mu)*(x-mu)) / (2*sigma*sigma));
	return y;
}
double dist_ParetoPDF (
	double x,
	double xm,
	double alpha)
{
	if (x < xm)
		return 0.0;
	return (alpha*pow (xm, alpha)) / pow (x, alpha+1);
}
double dist_ExponentialPDF (
	double x,
	double lambda)
{
	if (x < 0.0)
		return 0.0;
	return lambda * pow (NEPER, -lambda*x);
}

double dist_UnifCDF (
	double x,
	double a,
	double b)
{
	if (x < a)
		return a;
	else if (x > b)
		return b;
	return (x - a) / (b - a);
}
double dist_GaussCDF (
	double x,
	double mu,
	double sigma)
{
	return (1 + erf ((x - mu) / sqrt (2*sigma*sigma))) / 2;
}
double dist_ParetoCDF (
	double x,
	double xm,
	double alpha)
{
	if (x < 1.0)
		return 0.0;
	return 1.0 - (pow (xm, alpha) / pow (x, alpha));  
}
double dist_ExponentialCDF (
	double x,
	double lambda)
{
	if (x < 0.0)
		return 0.0;
	return 1.0 - pow (NEPER, -lambda*x); 
}


double dist_UnifInv (
	double a,
	double b)
{
	return (b - a) 
		* rng_MT19937_genrand () 
		+ a;
}

double dist_GaussBoxMuller ()
{
	static int f = OFF;
	static double u1, u2;
	f = f == ON ? OFF : ON;
	if (f == ON) {
		u1 = rng_MT19937_genrand ();
		u2 = rng_MT19937_genrand ();
	}
	return f == ON
		? sqrt (-2.0 * log (u1)) * sin (2.0 * PI * u2)
		: sqrt (-2.0 * log (u1)) * cos (2.0 * PI * u2);
}
double dist_GaussPolar ()
{
	static int f = OFF;
	static double u1, u2;
	static double w = DBL_MAX;
	f = f == ON ? OFF : ON;
	if (f == ON) {
		w = DBL_MAX;
		while (w >= 1.0 || w == 0.0) {
			u1 = 2 * rng_MT19937_genrand () - 1;
			u2 = 2 * rng_MT19937_genrand () - 1;
			w = u1*u1 + u2*u2;
		}
	}
	return f == ON
		? u1 * sqrt (-2 * log (w) / w)
		: u2 * sqrt (-2 * log (w) / w);
}

double dist_ParetoInv (
	double xm,
	double alpha)
{
	return xm / pow (1.0 - rng_MT19937_genrand (), 1.0 / alpha);
}

double dist_ExponentialInv (
	double lambda)
{
	return -log (1 - rng_MT19937_genrand ()) / lambda;
}


double* dist_CDFIncremental (
	proc_Points *pt)
{
	if (pt == NULL ||
		proc_CheckPoints (pt) == ERR)
		io_PrintErr (ERR, "invalid parameters in"
			" dist_CDF");
	int p;
	double *y = (double*) util_MemMalloc (pt->size * sizeof (double));
	for (p=0; p<pt->size; p++)
		y[p] = (double) (p + 1) / pt->size;
	return y;
}


void dist_Buckets (
	proc_Points *pt,
	int 		b,
	double 		*x,
	double 		*y,
	double 		*buck,
	int 		pdf)
{
	if (pt == NULL || b <= 0 || x == NULL || y == NULL ||
		proc_CheckPoints (pt) == ERR)
		io_PrintErr (ERR, "invalid parameters in"
			" dist_CDFBuckets");
	int p, j;
	double min, max;

	min = DBL_MAX;
	max = -DBL_MAX;
	for (p=0; p<pt->size; p++) {
		if (pt->points[p] < min)
			min = pt->points[p];
		if (pt->points[p] > max)
			max = pt->points[p];
	}
	
	*buck = (max - min) / b;
	x[0] = min + *buck / 2;

	for (j=1; j<b; j++)
		x[j] = x[j-1] + *buck;
	*buck /= 2;

	for (p=0; p<pt->size; p++)
		for (j=0; j<b; j++)
			if (pt->points[p] < x[j] + *buck) {
				y[j] += 1.0;
				if (pdf == ON)
					break;
			}
	if (pdf != ON)
		for (j=0; j<b; j++)
			y[j] /= pt->size;
}


double dist_F (
	double Rsq,
	int    pt)
{
	double F = Rsq / ((1.0 - Rsq) / (double) (pt - 2));
	if (TestHVerbosity > TestH_MEDIUM && TestHPrintPlain == OFF)
		fprintf (stdout, "   F: %.2lf / (%.2lf / %d)  \t= %.6lf\n\n", 
			Rsq, 1 - Rsq, pt - 2, F);
	return F;
}

