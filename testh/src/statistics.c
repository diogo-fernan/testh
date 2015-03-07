
#include "process.h"
#include "statistics.h"
#include "io.h"
#include "utilities.h"

#include <stdlib.h>
#include <math.h>
#include <time.h>


void stat_Mean (pts)
	proc_Points *pts;
{
	if (proc_CheckPoints (pts) == ERR)
		io_PrintErr (ERR, "invalid proc_Points structure in"
			" stat_Mean");
	int p;
	double sum;
	for (p=0, sum=0.0; p<pts->size; p++) {
		sum += pts->points[p];
		// fprintf (stdout, "%lf\n", pts->points[p]);
	}
	pts->mean = (double) sum / pts->size;
	// fprintf (stdout, "mean %d %lf %lf\n", pts->size, sum, pts->mean);
}

void stat_AbsoluteMomentsTime (pts, exponent)
	proc_Points *pts;
	int exponent;
{
	if (proc_CheckPoints (pts) == ERR || exponent <= 0)
		io_PrintErr (ERR, "invalid parameters in"
			" stat_AbsoluteMomentsTime");
	int p;
	long double sum;
	stat_Mean (pts);
	for (p=0, sum=0.0; p<pts->size; p++)
		sum += pow (pts->points[p] - pts->mean, exponent);
	pts->amt = (double) sum / (double) pts->size;
	pts->moment = exponent;
	// fprintf (stdout, "stdDev %lf %lf", pts->mean, pts->stdDev);
}

void stat_Variance (pts)
	proc_Points *pts;
{
	if (proc_CheckPoints (pts) == ERR)
		io_PrintErr (ERR, "invalid parameters in"
			" stat_Variance");
	stat_AbsoluteMomentsTime (pts, 2);
}

void stat_StdDeviation (pts)
	proc_Points *pts;
{
	if (proc_CheckPoints (pts) == ERR)
		io_PrintErr (ERR, "invalid proc_Points structure in"
			" stat_StdDeviation");
	stat_Variance (pts);
	pts->stdDev = (double) sqrt (pts->amt);
}

double stat_Autovariance (pts, k)
	proc_Points *pts;
	int k;
{
	if (proc_CheckPoints (pts) == ERR || k < 0)
		io_PrintErr (ERR, "invalid proc_Points structure in"
			" stat_Autovariance");

	int p;
	double sum;
	stat_Mean (pts);
	for (p=0, sum=0.0; p<pts->size-k; p++)
		sum += (pts->points[p] - pts->mean) * 
			(pts->points[p + k] - pts->mean);
	return sum / (pts->size - k);
}

// http://shadow.eas.gatech.edu/~jean/paleo/Meko_Autocorrelation.pdf page 2
void stat_AutocorrelationFunction (proc, kmin, kmax, h)
	proc_Process *proc;
	int kmin;
	int kmax;
	double h;
{
	if (proc == NULL || proc->points == NULL ||
		proc_CheckPoints (proc->points) == ERR || 
		proc->points->typeOfSignal != TestH_fGn ||
		kmin <= 0 || kmin <= 0 || kmin >= kmax)
		io_PrintErr (ERR, "invalid parameters in"
			" stat_AutocorrelationFunction");
	
	clock_t c;
	util_TimeIt (&c);

	int k;
	double r, gamma;
	proc_Points *pts = proc->points;
	stat_Mean (pts);
	if (TestHVerbosity > TestH_NONE && TestHPrintPlain == OFF) {
		proc_PrintProcess (proc);
		fprintf (stdout, "\n\t     Autocorrelation Function\n");
		fprintf (stdout, "\tlag k\t ac\t       exact fGn ac\n");
	}
	for (k=kmin; k<=kmax; k++) {
		r = stat_Autovariance (pts, k) / stat_Autovariance (pts, 0);
		gamma = (pow (k-1, 2*h) - 2*pow (k, 2*h) + pow (k+1, 2*h)) / 2;
		if (TestHVerbosity > TestH_NONE && TestHPrintPlain == OFF) {
			fprintf (stdout, " %10d\t%.6lf\t%.6lf\n", 
				k, r, gamma);
		}
	}
	fprintf (stdout, "\n");

	util_TimeWr (&c);
}
