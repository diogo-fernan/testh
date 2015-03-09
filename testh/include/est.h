#ifndef __TESTH_ESTIMATORS__
#define __TESTH_ESTIMATORS__

#include "proc.h"

enum TestHEstimator { 
	TestH_RS	= 0x1,
	TestH_VT 	= 0x2,
	TestH_AMT	= 0x3,
	TestH_EBP	= 0x4
};

typedef enum TestHEstimator est;


double est_RescaledRangeStatistics (proc_Process *pr);

double est_VarianceTime (proc_Process *pr);
void est_AbsoluteMomentsTime (proc_Process *pr, int mom);

void est_EmbeddedBranchingProcess (proc_Process *pr, int K);

#endif