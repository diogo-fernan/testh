
#include "process.h"
#include "batteries.h"
#include "io.h"
#include "generators.h"
#include "estimators.h"
#include "utilities.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// test all estimators for a given generator
void batt_Standard (
	gen 	generator,
	int 	n,
	int 	sim,
	double 	h1,
	double 	h2,
	double  inc,
	proc_ScalesConfig *conf)
{
	if (gen_CheckGen (generator) == ERR || n <= 0 || sim <= 0 ||
		io_CheckH (h1) == ERR || io_CheckH (h2) == ERR ||
		h1 >= h2 || inc > h2 - h1 ||
		proc_CheckScalesConfig (conf) == ERR) 
		io_PrintErr (ERR, "invalid parameters in"
			" batt_Standard");

	verb verbbak = TestHVerbosity;
	int printbak = TestHPrintPlain;

	TestHVerbosity = TestH_LOW;
	TestHPrintPlain = OFF;

	clock_t c;
	util_TimeIt (&c);

	if (generator == TestH_Gauss) {
		h1 = h2 = 0.5;
		inc = 0.0;
	}
	
	io_PrintSep ();
	TestHPrintSep = OFF;
	fprintf (stdout, "%sBATTERY%s: Standard\n",
		CGRAY_BLUE, CRESET);
	fprintf (stdout, " - %d sims, H: [%.2lf, %.2lf] += %.4lf\n",
		sim, h1, h2, inc);
	gen_PrintHeader (generator);
	fprintf (stdout, " - %d points\n", n);
	proc_PrintHeaderScales (conf);
	TestHVerbosity = TestH_NONE;
	fprintf (stdout, "\n   H\t RS\t\t\t VT\n");

if (generator == TestH_AR) {
	#define RENEWAL 10000
	#define xm 		(double) 1.0
} else if (generator == TestH_fBmSGA || generator == TestH_3SPG) {
	#define S 		16
}

#define EST 2

	proc_Process *proc = NULL;

	int i, j;
	double h, s = (double) sim;
	double *e, *em, *es;

	e  = (double*) util_MemCalloc (EST, sizeof (double));
	em = (double*) util_MemCalloc (EST, sizeof (double));
	es = (double*) util_MemCalloc (EST, sizeof (double));

	for (h=h1; h<h2+EPSILON; h+=inc) {
		for (j=0; j<EST; j++)
			em[j] =	es[j] = 0.0;

		for (i=0; i<sim; i++) {
			switch (generator) {
				case TestH_Gauss:
					proc = gen_Gaussian (n, TestH_fGn);
					break;
				case TestH_AR:
					proc = gen_AggRenewal (n, h, RENEWAL, xm, TestH_fGn);
					break;
				case TestH_fBmSGA:
					proc = gen_fBmSequentialGenerationAlgorithm (n, h, S, 
						TestH_fGn);
					break;
				case TestH_3SPG:
					proc = gen_SimpleSelfSimilarProcessGenerator (n, h, S, 
						TestH_fGn);
					break;
				case TestH_Hosk:
					proc = gen_Hosking (n, h, TestH_fGn);
					break;
				
				default:
					break;
			}

			proc_CreateScales (proc, conf);

			e[0] = est_RescaledRangeStatistics (proc);
			e[1] = est_VarianceTime (proc);

			for (j=0; j<EST; j++) {
				em[j] += e[j];
				es[j] += e[j] * e[j];
			}

			proc_DeleteProcess (proc);
		}

		fprintf (stdout, "  %.2lf", h);
		for (j=0; j<EST; j++) {
			em[j] /= s;
			es[j] /= s;
			fprintf (stdout, "\t%.4lf (%.2e)", 
				em[j], sqrt (es[j] - em[j] * em[j])); 

		}
		fprintf (stdout, "\n");
		if (generator == TestH_Gauss)
			break;
	}

	TestHVerbosity = verbbak;
	TestHPrintPlain = printbak;

	long int mem = sizeof (e) * EST * 3;
	util_MemWr (mem);
	util_MemFree (e);
	util_MemFree (em);
	util_MemFree (es);
	util_TimeWr (&c);
}

// Test all estimators for a given process
void batt_Estimators (
	proc_Process *proc)
{

}


