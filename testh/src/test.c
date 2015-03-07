
#include "process.h"
#include "generators.h"
#include "estimators.h"
#include "distributions.h"
#include "statistics.h"
#include "rng.h"
#include "io.h"
#include "utilities.h"
#include "batteries.h"

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

	proc_Process *proc1 = gen_fBmSequentialGenerationAlgorithm (100, 0.8, S, 
				TestH_fGn);
	// proc_Process *proc2 = gen_AggRenewal (100, 0.8, RENEWAL, xm, TestH_fGn);

	proc_PrintPoints (proc1->points);
	// proc_PrintPoints (proc2->points);

	proc_CreateScales (proc1, TestH_INC, 2, 8, 2);
	proc_PrintScales (proc1);

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

	proc_Process *proc = NULL;

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
			// proc = gen_SimpleSelfSimilarProcessGenerator (102400, j, S, TestH_fGn);
			// proc = gen_fBmSequentialGenerationAlgorithm (102400, j, S, TestH_fGn);
			// proc = gen_AggRenewal (100000, j, RENEWAL, xm, TestH_fGn);
			proc = gen_Hosking (10000, j, TestH_fGn);

			// proc_PrintPoints (proc->points);

			proc_CreateScales (proc, TestH_POW, 2, 10, 2);
			// proc_CreateScales (proc, TestH_INC, 2, 10, 2);
			// proc_PrintScales (proc);
			// proc_PrintPoints (&proc->scales->scales[0]);

			h1 = est_VarianceTime (proc);
			h2 = est_RescaledRangeStatistics (proc);
			// fprintf (stdout, "H1: %lf H2: %lf\n", H1 / 100.0, H2);
			H1 += h1;
			H2 += h2;
			H1_s += h1 * h1;
			H2_s += h2 * h2;

			proc_DeleteProcess (proc);
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

	proc_Process *proc = NULL;
	// proc_Process *proc2 = NULL;
	proc_ScalesConfig *conf = NULL;

	// proc = gen_ReadFile ("../files/PPRMI.csvTU1.0", "PPRMI.csvTU1.0", TestH_fGn);
	// proc = gen_ReadFile ("../files/BC-pAug89.txt.sed", "BC_Aug89", TestH_fGn);
	// proc = gen_Gaussian (N, TestH_fGn);
	// proc = gen_AggRenewal (N, H, RENEWAL, xm, TestH_fGn);
	// proc = gen_fBmSequentialGenerationAlgorithm (N, H, S, TestH_fGn);
	// proc = gen_SimpleSelfSimilarProcessGenerator (N, H, S, TestH_fGn);

	// proc2 = gen_ExternGen (100, "ExtGen", ExtGen, TestH_fGn);
	// proc_CreateScales (proc2, MIN_SCALE, MAX_SCALE, STEP);
	// proc_PrintProcess (proc2);
	// proc_PrintPoints (proc2->points);

	// proc_PrintPoints (proc->points);
	// proc_FractionalGaussianNoise (proc);
	// proc_Normalize (proc->points);
	// proc_PrintPoints (proc->points);
	// proc_FractionalBrownianMotion (proc);
	// proc_PrintPoints (proc->points);

	// conf = proc_CreateScalesConfig (TestH_INC, 4, 10, 2);
	conf = proc_CreateScalesConfig (TestH_POW, 7, 11, 2);

	// proc_CreateScales (proc, conf);
	// proc_CreateScalesTime (proc, TMIN_SCALE, TMAX_SCALE, TMULTIPLIER);
	
	// proc_PrintProcessStruct (proc);
	// proc_PrintProcess (proc);
	// proc_PrintPoints (&proc->scales->scales[0]);
	// proc_PrintScales (proc);

	// est_VarianceTime (proc);
	// est_RescaledRangeStatistics (proc);
	// est_AbsoluteMomentsTime (proc, 2);
	// est_AbsoluteMomentsTime (proc, 3);
	// est_EmbeddedBranchingProcess (proc, 5);

	batt_Standard (TestH_fBmSGA, 10e5, 2, 0.50, 0.99, 0.05, conf);

	// stat_AutocorrelationFunction (proc, 1, 40, H);

	// util_MemWr (proc_SizeOfProcess (proc));
	// proc_DeleteScales  (proc);
	proc_DeleteProcess (proc);

	io_PrintDone ();
	return EXIT_SUCCESS;
}
#endif


double ExtGen (void) {
	return rng_MT19937_genrand ();
}
