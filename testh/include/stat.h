#ifndef __TESTH_STATISTICS__
#define __TESTH_STATISTICS__

#include "proc.h"
#include "dist.h"

// http://www.danielmiessler.com/blog/the-craziest-thing-youll-ever-learn-about-pi
// 39 digits is enough
// But double cuts precision
#define PI 		(double) 3.141592653589793238462643383279502884197
#define NEPER 	(double) 2.71828

void stat_Mean (proc_Points *pt);
void stat_AbsoluteMomentsTime (proc_Points *pt,	int e);
void stat_StdDeviation (proc_Points *pt);
void stat_Variance (proc_Points *pt);
double stat_Min (proc_Points *pt);
double stat_Max (proc_Points *pt);
double stat_Autovariance (proc_Points *pt, int k);
double stat_Covariance (int k, double h);
void stat_Autocorrelation (proc_Process *pr, int kmin, int kmax, double h);

void stat_GOFKolmogorovSmirnov (proc_Points *pt, dist d);
void stat_GOFChiSquare (proc_Points *pt, dist d);

#endif