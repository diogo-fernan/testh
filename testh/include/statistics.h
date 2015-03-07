#ifndef __TESTH_STATISTICS__
#define __TESTH_STATISTICS__

#include "process.h"

// http://www.danielmiessler.com/blog/the-craziest-thing-youll-ever-learn-about-pi
// 39 digits is enough
// But double cuts precision
#define PI (double) 3.141592653589793238462643383279502884197


void stat_Mean (proc_Points *points);
void stat_AbsoluteMomentsTime (proc_Points *points,
	int exponent);
void stat_StdDeviation (proc_Points *points);
void stat_AutocorrelationFunction (proc_Process *proc, int kmin, int kmax, 
	double h);

#endif