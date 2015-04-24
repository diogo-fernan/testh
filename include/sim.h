#ifndef __TESTH_SIMULATORS__
#define __TESTH_SIMULATORS__

#include "proc.h"
#include "gen.h"


proc_Process* sim_NetworkTraffic (gen g, int n, double h, 
	int B, double L, int Imin, const char *pkts);


#endif