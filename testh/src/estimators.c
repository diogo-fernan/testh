
#include "io.h"
#include "utilities.h"
#include "estimators.h"
#include "regression.h"
#include "distributions.h"
#include "statistics.h"

#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <sys/param.h>

#define TESTH_MAX 256

static char TestHEstFile[TESTH_MAX] = {'\0'}; 


static void est_Init (
	est 	estimator,
	int 	res,
	char 	*pname,
	int 	moment)
{
	char e1[50], e2[50];

	switch (estimator) {
		case TestH_RS:
			strcpy (e1, "RescaledRangeStatistics");
			strcpy (e2, "RS");
			break;
		case TestH_VT:
			strcpy (e1, "VarianceTime");
			strcpy (e2, "VT");
			break;
		case TestH_AMT:
			sprintf (e1, "AbsoluteMomentsTime %d", moment);
			strcpy (e2, "AMT");
			break;
		case TestH_EBP:				
			strcpy (e1, "EmbeddedBranchingProcess");
			strcpy (e2, "EBP");
			break;

		default:
			strcpy (e1, "undefined");
			strcpy (e2, "undefined");
			break;
	}

	if (TestHEstWrToFile == ON) {
		memset (TestHEstFile, '\0', TESTH_MAX);
		strncpy (TestHEstFile, pname, strlen (pname));
		strncat (TestHEstFile, ".", 1);
		strncat (TestHEstFile, e2, strlen (e2));
		io_FileClean (TestHEstFile);
	}

	if (TestHVerbosity > TestH_NONE && TestHPrintPlain == OFF) {
		fprintf (stdout, "%sESTIMATOR%s: %s -- %s\n",
			CGRAY_BLUE, CRESET, e1, e2);
		if (TestHEstWrToFile == ON && TestHPrintPlain == OFF)
			fprintf (stdout, " - EstWrToFile: ON\n");
		fprintf (stdout, "\n");
		if (res == ON && TestHVerbosity > TestH_MEDIUM && TestHPrintPlain == OFF)
			fprintf (stdout, "  scale\t\t ln(scale)\t ln(%s)\n", e2);
	}
}

static void est_PrintResult (
	int 	scale,
	double 	est_log)
{
	if (TestHVerbosity > TestH_MEDIUM) {
		if (TestHPrintPlain == OFF)
			fprintf (stdout, " %5d\t\t",
				scale);
		fprintf (stdout, "%6.2lf\t\t%6.2lf\n",
			log (scale), est_log);
	}
	if (TestHEstWrToFile == ON)
		io_FileWr (TestHEstFile, "a", log (scale), est_log);
}

static void est_PrintH (
	est 	estimator,
	double 	h,
	double 	m,
	int 	moment)
{
	if (TestHVerbosity > TestH_NONE && TestHEstPrintH == ON) {
		if (TestHPrintPlain == ON)
			fprintf (stdout, "\t%lf\n", h);
		else {
			fprintf (stdout, "   d:\t\t\t\t= %lf\n", 2 - h);
			switch (estimator) {
				case TestH_RS:
					fprintf (stdout, "   H:\t\t\t\t= %.6lf\n", h); 
					break;
				case TestH_AMT:
					fprintf (stdout, "   H: 1.0 + %.2lf / %d\t\t= %.6lf\n", 
						m, moment, h);
					break;
				default:
					break;
			}
		}
	}
}


static double est_RescaledRangeStatisticsAux (
	proc_Points *fodp,  // original points, first order differences process 
	proc_Points *pts)   // scale points
{    
	if (proc_CheckPoints (pts) == ERR || proc_CheckPoints (fodp) == ERR)
		io_PrintErr (ERR, "invalid proc_Points structure(s) in"
			" est_RescaledRangeStatisticsAux");

	int i, j, l, c;
	long double ds, min, max;
	long double ssum, rs_tot;

	for (i=0, rs_tot=0.0; i<pts->size; i++) {
		max = -DBL_MAX; min = DBL_MAX;
		ssum = 0.0;
		c = 0;
		for (j=0; j<pts->scale; j++) {
			ds = 0.0;
			for (l=pts->scale*i; l<=pts->scale*i+j && l < fodp->size; l++) {
				ds += fodp->points[l] - pts->points[i];
				if (j == pts->scale - 1) {
					c++;
					ssum += fodp->points[l] * fodp->points[l];
				}
			}
			// if (j == pts->scale-1 && abs (ds) != 0.0)
			// 	fprintf (stdout, "\t%d %d %d %Lf\n", l, pts->scale, j, ds);
			// fprintf (stdout, "%lf\n", ssum);
				// ssum = 0.0 many times on the normalized process of
				// the aggregation of renewal processes
				// below NaN is obtained because sqrt (0.0)
			if (ds < min) min = ds;
			if (ds > max) max = ds;
		}
		if (min == DBL_MAX)
			min = 0.0;
		if (max == -DBL_MAX)
			max = 0.0;
		long double stdDev = sqrt (ssum / (double) c - 
			pts->points[i] * pts->points[i]);
			// stdDev = 0.0 ? 1.0 : stdDev;
		rs_tot += (max - min) / stdDev;
	}

	return rs_tot / pts->size;
}

double est_RescaledRangeStatistics (
	proc_Process *proc)
{
	if (proc_CheckProcess (proc) == ERR)
		io_PrintErr (ERR, "invalid proc_Process structure in"
			" est_RescaledRangeStatistics");
	
	clock_t c;
	util_TimeIt (&c);
	
	int s;
	double *y, *x;
	proc_Scales *scales;
	proc_Points *points;

	scales = proc->scales;
	x = (double*) util_MemMalloc (scales->conf->no * sizeof (double));
	y = (double*) util_MemMalloc (scales->conf->no * sizeof (double));

	proc_PrintHeader (proc);
	est_Init (TestH_RS, ON, proc->name, -1);

	for (s=0; s<scales->conf->no; s++) {
		points = &scales->scales[s];
		
		x[s] = log (points->scale);
		y[s] = log (est_RescaledRangeStatisticsAux (proc->points, points));

		est_PrintResult (points->scale, y[s]);
	}

	reg_Linear *l = reg_LeastSquareMeans (y, x, scales->conf->no, ON);
	double Rsquared = reg_CoefficientOfDetermination (y, x, scales->conf->no);
	dist_F (Rsquared, scales->conf->no);

	double H = l->m;
	est_PrintH (TestH_RS, H, -1, -1);

	long int mem = sizeof (x) * scales->conf->no +
				sizeof (y) * scales->conf->no;
	util_MemWr (mem);
	util_MemFree (l);
	util_TimeWr (&c);
	return H;
}


static double est_AbsoluteMomentsTimeAux (
	proc_Process 	*proc,
	int 			moment)
{
	clock_t c;
	util_TimeIt (&c);

	int s;
	double *x, *y;
	proc_Scales *scales; 

	scales = proc->scales;
	x = (double*) util_MemMalloc (scales->conf->no * sizeof (double));
	y = (double*) util_MemMalloc (scales->conf->no * sizeof (double));
	
	// stat_AbsoluteMomentsTime (proc->points, moment);
	// double y_o = log (proc->points->amt);

	for (s=0; s<scales->conf->no; s++) {
		stat_AbsoluteMomentsTime (&scales->scales[s], moment);
		y[s] = log (scales->scales[s].amt);
		x[s] = log (scales->scales[s].scale);

		est_PrintResult (scales->scales[s].scale, y[s]);
		// fprintf (stdout, " %lf\n", (y[s] - y_o) / (2 * x[s]) + 1);
	}

	reg_Linear *l = reg_LeastSquareMeans (y, x, scales->conf->no, ON);
	double Rsquared = reg_CoefficientOfDetermination (y, x, scales->conf->no);
	dist_F (Rsquared, scales->conf->no);

	// H = 1 + beta / moment
	double H = 1.0 + l->m / moment;
	est_PrintH (TestH_AMT, H, l->m, moment);

	long int mem = sizeof (x) * scales->conf->no +
				sizeof (y) * scales->conf->no;
	util_MemWr (mem);
	util_MemFree (l);
	util_TimeWr (&c);
	return H;
}


double est_VarianceTime (
	proc_Process *proc)
{
	if (proc_CheckProcess (proc) == ERR)
		io_PrintErr (ERR, "invalid proc_Process structure in"
			" est_VarianceTime");

	proc_PrintHeader (proc);
	est_Init (TestH_VT, ON, proc->name, -1);

	return est_AbsoluteMomentsTimeAux (proc, 2);
}

void est_AbsoluteMomentsTime (
	proc_Process 	*proc,
	int 			moment)
{
	if (proc_CheckProcess (proc) == ERR || moment <= 0)
		io_PrintErr (ERR, "invalid parameters in"
			" est_AbsoluteMomentsTime");

	proc_PrintHeader (proc);
	est_Init (TestH_AMT, ON, proc->name, moment);

	est_AbsoluteMomentsTimeAux (proc, moment);
}


static double est_CrossingLines (
	double 	delta,
	int 	k)
{
	return (double) delta * pow (2, k);
}

void est_EmbeddedBranchingProcess (
	proc_Process 	*proc,
	int 			K)
{
	if (proc->points == NULL || proc_CheckPoints (proc->points) == ERR)
		io_PrintErr (ERR, "invalid parameters in"
			" est_EmbeddedBranchingProcess");
	
	clock_t c;
	util_TimeIt (&c);
	proc_PrintHeader (proc);
	est_Init (TestH_EBP, OFF, proc->name, -1);
 
	int p, k;
	int *crossings;
	double line, lastLineCrossed, delta;
	double *crossingLines;

	proc_Points *fGmPoints = proc->points;

	if (proc->points->typeOfSignal == TestH_fGn) {
		delta = proc->points->stdDev;
		// fGmPoints = proc_fractionalBrownianMotion (proc->points);
		// not supposed to be commented
	}
	else if (proc->points->typeOfSignal == TestH_fBm) {
		proc_Points* fGnPoints = NULL; // = proc_fractionalGaussianNoise (proc->points);
		// not supposed to be commented
		delta = fGnPoints->stdDev;
		util_MemFree (fGnPoints);
	}

	stat_StdDeviation (fGmPoints);

	lastLineCrossed = DBL_MAX;
	crossings = (int*) util_MemCalloc (K, sizeof (int));
	crossingLines = (double*) util_MemMalloc (K * sizeof (double));

	for (k=1; k<=K; k++) {
		crossingLines[k-1] = est_CrossingLines(2, k);
		fprintf (stdout, "line %d: %lf\n", k, crossingLines[k-1]);
	}

	fprintf (stdout, "\n");
	for (p=1; p<fGmPoints->size; p++) {
		for (k=1; k<=K; k++) {
			line = crossingLines[k-1];
			if (((fGmPoints->points[p-1] <= line &&
					fGmPoints->points[p] >= line) ||
				(fGmPoints->points[p-1] >= line &&
					fGmPoints->points[p] <= line)) &&
				line != lastLineCrossed) {
					
					fprintf (stdout, "point %d crossed line %lf: %lf --> %lf\n",
						p, line, fGmPoints->points[p-1], fGmPoints->points[p]);
	
				lastLineCrossed = line;
				crossings[k]++;
				break;
			}
		}
	}

	fprintf (stdout, "\n");
	int cardinalitySum = 0;
	double uk, u = 0.0;
	double *x, *y;
	
	x = (double*) util_MemMalloc (K * sizeof (double));
	y = (double*) util_MemMalloc (K * sizeof (double));

	for (k=0; k<K; k++) {
		if (k < K - 1)
			uk = (double) crossings[k] / crossings[k+1];
		else
			uk = (double) crossings[k] / crossings[k];
		fprintf (stdout, "u_%d: %lf = %d / %d\n",
			k+1, uk, crossings[k], crossings[k+1]);
		
		u += (double) crossings[k] * uk;
		cardinalitySum += crossings[0];
		
		x[k] = crossings[k] / crossings[K-1];
		y[k] = pow (2, k + 1);
	}

	u /= cardinalitySum;
	fprintf (stdout, "\n H = log (2) / log (%lf) = %lf\n",
		u, log (2) / log (u));

	reg_LeastSquareMeans (y, x, K, ON);
	double Rsquared = reg_CoefficientOfDetermination (y, x, K);
	dist_F (Rsquared, K);
	
	util_TimeWr (&c);
}

