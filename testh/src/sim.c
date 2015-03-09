
#include "proc.h"
#include "gen.h"
#include "dist.h"
#include "sim.h"
#include "stat.h"
#include "util.h"
#include "rng.h"
#include "io.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>





proc_Process* sim_NetworkTraffic (
	gen 	g,
	int 	n,
	double 	h,
	int 	B,
	double 	L,
	int 	Imin,
	const char *f)
{
	if (gen_CheckGen (g) == ERR || n <= 0 || 
		io_CheckH (h) == ERR || B <= 0 || L <= 0.0 || L > 1.0 || Imin <= 0) 
		io_PrintErr (ERR, "invalid parameters in"
			" sim_NetworkTraffic");

	verb verbbak = TestHVerbosity;
	int  bak1 = TestHPrintPlain;
	int  bak2 = TestHPrintSep;

	TestHVerbosity = TestH_LOW;
	TestHPrintPlain = OFF;

	clock_t c;
	util_TimeIt (&c);

	io_PrintSep ();
	TestHPrintSep = OFF;
	fprintf (stdout, "%sSIMULATOR%s: NetworkTraffic\n",
		CGRAY_BLUE, CRESET);
	gen_PrintHeader (g);
	fprintf (stdout, " - %d points\n", n);
	TestHVerbosity = TestH_NONE;

	int p, j;
	double r, np, Bu, Bnu, Iavg, Istddev, I, P;
	proc_Process *pr = NULL;

	proc_Process *pkts = gen_ReadFile (f, NULL, TestH_fBm);
	double *y = dist_CDFIncremental (pkts->points);
	pr = gen_GenProc (g, n, h);

	stat_Mean (pkts->points);
	Bu 	= L * B;
	Bnu	= B - Bu;
	np 	= Bu / pkts->points->mean;
	Iavg = Bnu / np;
	Istddev = (Iavg - Imin) / 2;

	pr->points->times = (double*) util_MemMalloc (n * sizeof (double));

	for (p=0; p<n; p++) {
		r = rng_MT19937_genrand ();
		for (j=0; j<pkts->points->size; j++)
			if (r < y[j])
				break;

		P = pkts->points->points[j];
		I = pr->points->points[p] * Istddev + Iavg;
		if (I < Imin)
			I = Imin;
		else if (I > 2 * Iavg - Imin)
			I = 2 * Iavg - Imin;
		I /= Bu;
		// I *= 10;

		pr->points->points[p] = P;
		if (p > 0)
			pr->points->times[p] = pr->points->times[p-1] + I;
		else
			pr->points->times[p] = I;
	}

	TestHVerbosity = verbbak;
	TestHPrintPlain = bak1;
	TestHPrintSep 	= bak2;

	long int mem = 0.0;
	util_MemWr (mem);
	// util_MemFree (. . .);
	util_TimeWr (&c);
	return pr;
}


