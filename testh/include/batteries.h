#ifndef __TESTH_BATTERIES__
#define __TESTH_BATTERIES__

#include "process.h"
#include "generators.h"

#define EPSILON 0.000000000001


void batt_Standard (gen generator, int n, int sim, double h1, double h2, 
	double inc, proc_ScalesConfig *conf);

#endif