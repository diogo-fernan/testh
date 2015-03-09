
#include "proc.h"
#include "gen.h"
#include "est.h"
#include "dist.h"
#include "stat.h"
#include "rng.h"
#include "io.h"
#include "util.h"
#include "batt.h"
#include "sim.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define N			100000
#define H 			(double) 0.8
#define S           16

#define RENEWAL		10000
#define xm 			(double) 1.0

#define MIN_SCALE	2
#define MAX_SCALE	20
#define STEP		2
#define N_SCALES	((MAX_SCALE - MIN_SCALE) / STEP + 1)

#define TMIN_SCALE  (double) 0.001
#define TMAX_SCALE  10
#define TMULTIPLIER 10


double ExtGen ();


#if 0

	// Another main here
	// set ifs to 0 on the others or change function name (e.g., to main2)

int main ()
{
	TestHVerbosity = TestH_HIGH;
	TestHEstWrToFile = ON;

	int i;
	double a[1000000];
	for (i=0; i<1000000; i++)
		a[i] = dist_Normal ();
		// fprintf (stdout, "%lf\n", dist_Normal ());
	proc_Process *p = proc_CreateProcess ("", a, NULL, 1000000, TestH_fGn);
	proc_CreateScales (p, TestH_INC, 10, 20, 2);
	est_VarianceTime (p);
	est_RescaledRangeStatistics (p);
}



#endif


#if 0
int main (int argc, char *argv[]) {
	TestHVerbosity = TestH_NONE;
	TestHPrintSep = OFF;
	TestHPrintHeader = OFF;
	TestHPrintMemCPU = OFF;
	TestHEstWrToFile = OFF;

	io_PrintInit ("make; make run", argv[0]);

	proc_Process *pr1 = gen_fBmSequentialGenerationAlgorithm (100, 0.8, S, 
				TestH_fGn);
	// proc_Process *pr2 = gen_AggRenewal (100, 0.8, RENEWAL, xm, TestH_fGn);

	proc_PrintPoints (pr1->points);
	// proc_PrintPoints (pr2->points);

	proc_CreateScales (pr1, TestH_INC, 2, 8, 2);
	proc_PrintScales (pr1);

	io_PrintDone ();
	return EXIT_SUCCESS;
}
#endif


#if 0
int main (int argc, char *argv[]) {
	TestHVerbosity = TestH_NONE;
	TestHPrintSep = OFF;
	TestHPrintHeader = OFF;
	TestHPrintMemCPU = OFF;
	TestHEstWrToFile = OFF;

	io_PrintInit ("make; make run", argv[0]);

	proc_Process *pr = NULL;

	int i;
	double j;
	double s = 10.0;
	double h1, h2, H1, H2;
	double H1_s, H2_s;

	fprintf (stdout, "\n  H & VT & R/S\n");
	for (j=0.50; j<1.0; j+=0.05) {
		H1 = H2 = 0.0;
		H1_s = H2_s = 0.0;
		for (i=0; i<(int) s; i++) {
			// pr = gen_SimpleSelfSimilarProcessGenerator (102400, j, S, TestH_fGn);
			// pr = gen_fBmSequentialGenerationAlgorithm (102400, j, S, TestH_fGn);
			// pr = gen_AggRenewal (100000, j, RENEWAL, xm, TestH_fGn);
			pr = gen_Hosking (10000, j, TestH_fGn);

			// proc_PrintPoints (pr->points);

			proc_CreateScales (pr, TestH_POW, 2, 10, 2);
			// proc_CreateScales (pr, TestH_INC, 2, 10, 2);
			// proc_PrintScales (pr);
			// proc_PrintPoints (&pr->scales->scales[0]);

			h1 = est_VarianceTime (pr);
			h2 = est_RescaledRangeStatistics (pr);
			// fprintf (stdout, "H1: %lf H2: %lf\n", H1 / 100.0, H2);
			H1 += h1;
			H2 += h2;
			H1_s += h1 * h1;
			H2_s += h2 * h2;

			proc_DeleteProcess (pr);
		}
		H1 /= s;
		H1_s /= s;
		H2 /= s;
		H2_s /= s;
		fprintf (stdout, " %.2lf & %lf (%lf) & %lf (%lf)\n", 
			j, H1, sqrt (H1_s - H1 * H1), H2, sqrt (H2_s - H2 * H2));
	}

	io_PrintDone ();
	return EXIT_SUCCESS;
}
#endif

#if 1
int main (int argc, char *argv[]) {
	TestHVerbosity = TestH_HIGH;
	TestHPrintPlain = OFF;
	TestHPrintSep = ON;
	TestHPrintHeader = ON;
	TestHPrintMemCPU = ON;
	TestHEstPrintH = ON;
	TestHEstWrToFile = OFF;

	io_PrintInit ("make; make run", argv[0]);
	// io_PrintTestH ();

	proc_Process *pr = NULL;
	// proc_Process *pr2 = NULL;
	proc_ScalesConfig *conf = NULL;

	// pr = gen_ReadFile ("../files/PPRMI.csvTU1.0", "PPRMI.csvTU1.0", TestH_fGn);
	// pr = gen_ReadFile ("../files/BC-pAug89.txt.sed", "BC_Aug89", TestH_fGn);
	// pr = gen_ReadFileMatrix ("data", "name", TestH_fGn);
	// pr = gen_Gaussian (N, TestH_fGn);
	// pr = gen_AggRenewal (N, H, RENEWAL, xm, TestH_fGn);
	// pr = gen_fBmSequentialGenerationAlgorithm (N, H, S, TestH_fGn);
	// pr = gen_SimpleSelfSimilarProcessGenerator (N, H, S, TestH_fGn);

	// pr2 = gen_ExternGen (100, "ExtGen", ExtGen, TestH_fGn);
	// proc_CreateScales (pr2, MIN_SCALE, MAX_SCALE, STEP);
	// proc_PrintProcess (pr2);
	// proc_PrintPoints (pr2->points);

	// proc_PrintPoints (pr->points);
	// proc_FractionalGaussianNoise (pr);
	// proc_Normalize (pr->points);
	// proc_PrintPoints (pr->points);
	// proc_FractionalBrownianMotion (pr);
	// proc_PrintPoints (pr->points);

	// conf = proc_CreateScalesConfig (TestH_INC, 4, 10, 2);
	conf = proc_CreateScalesConfig (TestH_POW, 7, 11, 2);

	// proc_CreateScales (pr, conf);
	// proc_CreateScalesTime (pr, TMIN_SCALE, TMAX_SCALE, TMULTIPLIER);
	
	// proc_PrintProcessStruct (pr);
	// proc_PrintProcess (pr);
	// proc_PrintPoints (&pr->scales->scales[0]);
	// proc_PrintScales (pr);

	// est_VarianceTime (pr);
	// est_RescaledRangeStatistics (pr);
	
	batt_Standard (TestH_fBmSGA, 10e5, 2, 0.50, 0.99, 0.05, conf);

	// stat_AutocorrelationFunction (pr, 1, 40, H);

	// util_MemWr (proc_SizeOfProcess (pr));
	// proc_DeleteScales (pr);
	proc_DeleteProcess (pr);

	io_PrintDone ();
	return EXIT_SUCCESS;
}
#endif


double ExtGen (void) {
	return rng_MT19937_genrand ();
}
