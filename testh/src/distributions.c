
#include "statistics.h"
#include "distributions.h"
#include "rng.h"
#include "io.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

double dist_Unif (a, b)
	double a;
	double b;
{
	return (b - a) 
		* rng_MT19937_genrand () 
		+ a;
}

double dist_NormalBoxMuller ()
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
double dist_NormalPolar ()
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

double dist_ParetoInv (xm, h)
	double xm;
	double h;
{
	return xm / pow(1.0 - rng_MT19937_genrand (), 1.0 /
		(3.0 - 2.0 * h)); // alpha
}

double dist_F (Rsquared, points)
	double Rsquared;
	int points;
{
	double F = Rsquared / ((1.0 - Rsquared) / (double) (points - 2));
	if (TestHVerbosity > TestH_MEDIUM && TestHPrintPlain == OFF)
		fprintf (stdout, "   F: %.2lf / (%.2lf / %d)  \t= %.6lf\n\n", 
			Rsquared, 1 - Rsquared, points - 2, F);
	return F;
}










