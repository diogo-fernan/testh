
#include "io.h"
#include "utilities.h"
#include "regression.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>



double reg_LinearModel (
	reg_Linear 	*l,
	double 		x)
{
	if (l == NULL)
		io_PrintErr (ERR, "invalid reg_Linear pointer (NULL) in"
			" reg_LinearModel");
	// y = mx + b
	return l->m * x + l->b;
}



double reg_CoefficientOfDetermination (
	double 	*y,
	double 	*x,
	int 	n)
{
	if (n <= 0 || x == NULL || y == NULL)
		io_PrintErr (ERR, "invalid parameters in"
			" reg_CoefficientOfDetermination");

	int i;
	double aux, y_mean, y_tot, SSres, SStot;

	y_tot = SSres = SStot = 0.0;
	reg_Linear *linear = reg_LeastSquareMeans (y, x, n, OFF);

	for (i=0; i<n; i++) {
		y_tot += y[i]; 
	}
	y_mean = y_tot / n;
	for (i=0; i<n; i++) {
		SStot += (y[i] - y_mean) * (y[i] - y_mean);
		aux = y[i] - reg_LinearModel (linear, x[i]);
		SSres += aux * aux;
	}
	double Rsquared = 1.0 - (SSres / SStot);
	if (TestHVerbosity > TestH_MEDIUM && TestHPrintPlain == OFF)
		fprintf (stdout, " R^2: 1.0 - (%.2lf / %.2lf)\t= %.6lf\n",
			SSres, SStot, Rsquared);
	return Rsquared;
}

// http://hotmath.com/hotmath_help/topics/line-of-best-fit.html
reg_Linear* reg_LeastSquareMeans (
	double 	*y,
	double 	*x,
	int 	n,
	int 	print)
{
	if (n <= 0 || x == NULL || y == NULL)
		io_PrintErr (ERR, "invalid parameters in"
			" reg_LeastSquareMethod");
	
	int i;
	double sx, sxx, sy, sxy;
	
	sx = sxx = sy = sxy = 0.0;

	for (i=0; i<n; i++) {
		sx  += x[i];	
		sxx += x[i] * x[i];
		sy  += y[i];
		sxy += x[i] * y[i];
	}

	sx 	/= (double) n;
	sxx /= (double) n;
	sy 	/= (double) n;
	sxy	/= (double) n;

	reg_Linear *r = (reg_Linear*) util_MemMalloc (sizeof (reg_Linear));

	r->m = (sxy - sx * sy) / (sxx - sx * sx);
	r->b = sy - r->m * sx;

	if (TestHVerbosity > TestH_MEDIUM && TestHPrintPlain == OFF && print == ON)
		fprintf (stdout, "\n LSM: y = %.2lf * x + %.2lf\n", r->m, r->b);
	return r;
}



