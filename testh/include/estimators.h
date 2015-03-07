#ifndef __TESTH_ESTIMATORS__
#define __TESTH_ESTIMATORS__

#include "process.h"

enum TestHEstimator { 
	TestH_RS	= 0x1,
	TestH_VT 	= 0x2,
	TestH_AMT	= 0x3,
	TestH_EBP	= 0x4
};

typedef enum TestHEstimator est;


double est_RescaledRangeStatistics (proc_Process *proc);

double est_VarianceTime (proc_Process *proc);
void est_AbsoluteMomentsTime (proc_Process *proc, int moment);

void est_EmbeddedBranchingProcess (proc_Process *proc, int K);

#endif