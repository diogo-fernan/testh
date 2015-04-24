
#include "proc.h"
#include "stat.h"
#include "io.h"
#include "util.h"
#include "dist.h"
#include "func.h"

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>


void stat_Mean (
	proc_Points *pt)
{
	if (proc_CheckPoints (pt) == ERR)
		io_PrintErr (ERR, "invalid proc_Points structure in"
			" stat_Mean");
	int p;
	double sum;
	for (p=0, sum=0.0; p<pt->size; p++)
		sum += pt->points[p];
	pt->mean = (double) sum / pt->size;
}

void stat_AbsoluteMomentsTime (
	proc_Points *pt,
	int 		e)
{
	if (proc_CheckPoints (pt) == ERR || e <= 0)
		io_PrintErr (ERR, "invalid parameters in"
			" stat_AbsoluteMomentsTime");
	int p;
	long double sum;
	stat_Mean (pt);
	for (p=0, sum=0.0; p<pt->size; p++)
		sum += pow (pt->points[p] - pt->mean, e);
	pt->amt = (double) sum / (double) pt->size;
	pt->moment = e;
}

void stat_Variance (
	proc_Points *pt)
{
	if (proc_CheckPoints (pt) == ERR)
		io_PrintErr (ERR, "invalid parameters in"
			" stat_Variance");
	stat_AbsoluteMomentsTime (pt, 2);
}

void stat_StdDeviation (
	proc_Points *pt)
{
	if (proc_CheckPoints (pt) == ERR)
		io_PrintErr (ERR, "invalid proc_Points structure in"
			" stat_StdDeviation");
	stat_Variance (pt);
	pt->stdDev = (double) sqrt (pt->amt);
}

double stat_Min (
	proc_Points *pt)
{
	if (proc_CheckPoints (pt) == ERR)
		io_PrintErr (ERR, "invalid proc_Points structure in"
			" stat_Min");
	int p, min = pt->points[0];
	for (p=1; p<pt->size; p++)
		if (pt->points[p] < min)
			min = pt->points[p];
	return min;
}
double stat_Max (
	proc_Points *pt)
{
	if (proc_CheckPoints (pt) == ERR)
		io_PrintErr (ERR, "invalid proc_Points structure in"
			" stat_Max");
	int p, max = pt->points[0];
	for (p=1; p<pt->size; p++)
		if (pt->points[p] > max)
			max = pt->points[p];
	return max;
}

double stat_Autovariance (
	proc_Points *pt,
	int 		k)
{
	if (proc_CheckPoints (pt) == ERR || k < 0)
		io_PrintErr (ERR, "invalid proc_Points structure in"
			" stat_Autovariance");
	int p;
	double sum;
	stat_Mean (pt);
	for (p=0, sum=0.0; p<pt->size-k; p++)
		sum += (pt->points[p] - pt->mean) * 
			(pt->points[p + k] - pt->mean);
	return sum; // / (pt->size - k);
}

double stat_Covariance (
	int 	k,
	double 	h)
{
	if (k == 0)
		return 1.0;
	else
		return (pow (k-1, 2*h) - 2 * pow (k, 2*h) + pow (k+1, 2*h)) / 2;
}

// http://shadow.eas.gatech.edu/~jean/paleo/Meko_Autocorrelation.pdf page 2
void stat_Autocorrelation (
	proc_Process 	*pr,
	int 			kmin,
	int 			kmax,
	double 			h)
{
	if (pr == NULL || pr->points == NULL ||
		proc_CheckPoints (pr->points) == ERR || 
		pr->points->typeOfSignal != TestH_fGn ||
		kmin <= 0 || kmin <= 0 || kmin >= kmax)
		io_PrintErr (ERR, "invalid parameters in"
			" stat_Autocorrelation");
	
	clock_t c;
	util_TimeIt (&c);

	int k;
	double r, gamma;
	proc_Points *pts = pr->points;
	stat_Mean (pts);
	if (TestHVerbosity > TestH_NONE && TestHPrintPlain == OFF) {
		proc_PrintProcess (pr);
		fprintf (stdout, "\n\t     Autocorrelation Function\n");
		fprintf (stdout, "\tlag k\t ac\t       exact fGn ac\n");
	}
	for (k=kmin; k<=kmax; k++) {
		r = stat_Autovariance (pts, k) / stat_Autovariance (pts, 0);
		gamma = stat_Covariance (k, h);
		if (TestHVerbosity > TestH_NONE && TestHPrintPlain == OFF) {
			fprintf (stdout, " %10d\t%.6lf\t%.6lf\n", 
				k, r, gamma);
		}
	}
	fprintf (stdout, "\n");

	util_TimeWr (&c);
}

#define SIGNIF_LEVEL	0.05
#define BINS1 			10
#define BINS2 			20
#define KS_MAX_K 		100

void stat_GOFKolmogorovSmirnov (
	proc_Points *pt, 
	dist 		d)
{
	if (pt == NULL || proc_CheckPoints (pt) == ERR ||
		dist_CheckDist (d) == ERR)
		io_PrintErr (ERR, "invalid parameters in"
			" stat_GOFKolmogorovSmirnov");

	int p;
	double prob, sup = -DBL_MAX;
	double *x, *y, buck;

	if (TestHVerbosity > TestH_NONE && TestHPrintPlain == OFF) {
		io_PrintSep ();
		fprintf (stdout, "GOF TEST: Kolmogorov-Smirnov");
		fprintf (stdout, "\n - %llu points\n", pt->size);
		fprintf (stdout, "\n\tsignif level:\t%.3lf\n\tbins:\t\t%d\n"
			"\tKS max:\t\t%d\n", SIGNIF_LEVEL, BINS2, KS_MAX_K);
	}

	x = (double*) util_MemMalloc (BINS2 * sizeof (double));
	y = (double*) util_MemCalloc (BINS2, sizeof (double));
	dist_Buckets (pt, BINS2, x, y, &buck, OFF);
	// x = pt->points;
	// y = dist_CDFIncremental (pt);

	for (p=0; p<BINS2; p++) {
		prob = dist_CDFProb (d, x[p]);
		if (fabs (y[p] - prob) > sup)
			sup = fabs (y[p] - prob);
	}
	sup *= sqrt (BINS2);

	int k;
	double P, sum;
	double term = -2 * sup*sup;

	for (k=1, sum=0.0; k<=KS_MAX_K; k++)
		// KS CDF
		sum += pow (-1, k-1) * exp (k*k * term);
	P = 1 - 2*sum;

	if (TestHVerbosity > TestH_NONE && TestHPrintPlain == OFF) {
		fprintf (stdout, "\n\tp-value:\t%lf\n", P);
		if (P > SIGNIF_LEVEL)
			fprintf (stdout, "\t\treject\n");
		else
			fprintf (stdout, "\t\taccept\n");
	}
}

void stat_GOFChiSquare (
	proc_Points *pt, 
	dist 		d)
{
	if (pt == NULL || proc_CheckPoints (pt) == ERR ||
		dist_CheckDist (d) == ERR)
		io_PrintErr (ERR, "invalid parameters in"
			" stat_GOFChiSquare");

	int p;
	double chsq, pr, P;
	double *x, *y, buck;

	if (TestHVerbosity > TestH_NONE && TestHPrintPlain == OFF) {
		io_PrintSep ();
		fprintf (stdout, "GOF TEST: Chi-Square");
		fprintf (stdout, "\n - %llu points\n", pt->size);
		fprintf (stdout, "\n\tsignif level:\t%.3lf\n\tbins:\t\t%d\n"
			"\tdf:\t\t%d\n", SIGNIF_LEVEL, BINS1, BINS1-1);
	}

	x = (double*) util_MemMalloc (BINS1 * sizeof (double));
	y = (double*) util_MemCalloc (BINS1, sizeof (double));
	dist_Buckets (pt, BINS1, x, y, &buck, ON);

	for (p=0, chsq=0.0; p<BINS1; p++) {
		pr = func_qsimp (d, x[p]-buck, x[p]+buck);
		if (pr != 0.0) {
			pr *= pt->size;
			chsq += ((y[p] - pr) * (y[p] - pr)) / pr;
		}
	}

	P = func_gammq (0.5 * (BINS1 - 1), 0.5 * chsq);
	fprintf (stdout, "\n\tp-value:\t%.30lf\n", P);
}

