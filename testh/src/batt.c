
#include "proc.h"
#include "batt.h"
#include "io.h"
#include "gen.h"
#include "est.h"
#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define EPSILON (double) 1.0e-12
#define EST 	2


void batt_Generator (
	gen 	g,
	int 	n,
	int 	sim,
	double 	h1,
	double 	h2,
	double  inc,
	proc_ScalesConfig *conf)
{
	if (gen_CheckGen (g) == ERR || n <= 0 || sim <= 0 ||
		io_CheckH (h1) == ERR || io_CheckH (h2) == ERR ||
		h1 >= h2 || inc > h2 - h1 ||
		proc_CheckScalesConfig (conf) == ERR) 
		io_PrintErr (ERR, "invalid parameters in"
			" batt_Generator");

	verb verbbak = TestHVerbosity;
	int  bak1 = TestHPrintPlain;
	int  bak2 = TestHPrintSep;

	TestHVerbosity = TestH_LOW;
	TestHPrintPlain = OFF;

	clock_t c;
	util_TimeIt (&c);

	if (g == TestH_Gaussian) {
		h1 = h2 = 0.5;
		inc = 0.0;
	}
	
	io_PrintSep ();
	TestHPrintSep = OFF;
	fprintf (stdout, "%sBATTERY%s: Generator\n",
		CGRAY_BLUE, CRESET);
	fprintf (stdout, " - %d sims, H: [%.2lf, %.2lf] += %.4lf\n",
		sim, h1, h2, inc);
	gen_PrintHeader (g);
	fprintf (stdout, " - %d points\n", n);
	proc_PrintHeaderScales (conf);
	TestHVerbosity = TestH_NONE;
	fprintf (stdout, "\n   H\t RS\t\t\t VT\n");

	proc_Process *pr = NULL;

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
			pr = gen_GenProc (g, n, h);
			proc_CreateScales (pr, conf);

			e[0] = est_RescaledRangeStatistics (pr);
			e[1] = est_VarianceTime (pr);

			for (j=0; j<EST; j++) {
				em[j] += e[j];
				es[j] += e[j] * e[j];
			}

			proc_DeleteProcess (pr);
		}

		fprintf (stdout, "  %.2lf", h);
		for (j=0; j<EST; j++) {
			em[j] /= s;
			es[j] /= s;
			fprintf (stdout, "\t%.4lf (%.2e)", 
				em[j], sqrt (es[j] - em[j] * em[j])); 

		}
		fprintf (stdout, "\n");
		if (g == TestH_Gaussian)
			break;
	}

	TestHVerbosity 	= verbbak;
	TestHPrintPlain = bak1;
	TestHPrintSep 	= bak2;

	long int mem = sizeof (e) * EST * 3;
	util_MemWr (mem);
	util_MemFree (e);
	util_MemFree (em);
	util_MemFree (es);
	util_TimeWr (&c);
}


void batt_Process (
	proc_Process 		*pr,
	proc_ScalesConfig 	*conf)
{
	if (pr == NULL ||
		proc_CheckPoints (pr->points) == ERR ||
		proc_CheckScalesConfig (conf) == ERR) 
		io_PrintErr (ERR, "invalid parameters in"
			" batt_Process");

	clock_t c;
	util_TimeIt (&c);

	double *e = (double*) util_MemCalloc (EST, sizeof (double));
	verb bak1 = TestHVerbosity;
	int  bak2 = TestHPrintHeader;
	int  bak3 = TestHPrintSep;
	int  bak4 = TestHPrintMemCPU;
	int  bak5 = TestHPrintPlain;

	TestHVerbosity 	 = TestH_HIGH;
	TestHPrintPlain  = OFF;
	TestHPrintMemCPU = OFF;
	io_PrintSep ();
	TestHPrintSep = OFF;
	fprintf (stdout, "%sBATTERY%s: Process\n", CGRAY_BLUE, CRESET);
	fprintf (stdout, " - process: %s\n", pr->name);
	fprintf (stdout, " - %llu points\n", pr->points->size);
	proc_PrintHeaderScales (conf);
	TestHVerbosity = TestH_NONE;

	proc_CreateScales (pr, conf);

	TestHPrintHeader = OFF;
	TestHVerbosity = bak1;

	e[0] = est_RescaledRangeStatistics (pr);
	fprintf (stdout, "\n");
	e[1] = est_VarianceTime (pr);

	fprintf (stdout, "\n\t  Results\n\tEst\t H\n");
	fprintf (stdout, "\t RS\t%.4lf\n", e[0]); 
	fprintf (stdout, "\t VT\t%.4lf\n", e[1]); 

	TestHPrintHeader = bak2;
	TestHPrintSep 	 = bak3;
	TestHPrintMemCPU = bak4;
	TestHPrintPlain  = bak5;

	long int mem = sizeof (e) * EST;
	util_MemWr (mem);
	util_MemFree (e);
	util_TimeWr (&c);
}


