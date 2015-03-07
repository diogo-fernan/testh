#ifndef __TESTH_GENERATORS__
#define __TESTH_GENERATORS__

#include "process.h"

enum TestHGenerator { 
	TestH_EG		= -0x1,
	TestH_RF		= -0x2,
	TestH_RFT		= -0x3,

	TestH_Gauss		= 0x0,

	TestH_AR 		= 0x1,
	TestH_fBmSGA	= 0x2,
	TestH_3SPG		= 0x3,
	TestH_Hosk		= 0x4
};

typedef enum TestHGenerator gen;


void gen_PrintHeader (gen generator);
int gen_CheckGen (gen generator);

proc_Process* gen_ExternGen (int N, char *name, double (*gen_func) (void), 
	tosig typeOfSignal);


proc_Process* gen_ReadFile (const char *filename, const char *name,
	tosig typeOfSignal);
proc_Process* gen_ReadFileTime (const char *filename, const char *name,
	tosig typeOfSignal);


proc_Process* gen_Gaussian (int N, tosig typeOfSignal);


proc_Process* gen_AggRenewal (int N, double h, int ren, int xm, 
	tosig typeOfSignal);


proc_Process* gen_fBmSequentialGenerationAlgorithm (int N, double h, 
	int scale_no, tosig typeOfSignal);
proc_Process* gen_SimpleSelfSimilarProcessGenerator ();


proc_Process* gen_Hosking (int N, double h, tosig typeOfSignal);

#endif