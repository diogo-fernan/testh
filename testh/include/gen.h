#ifndef __TESTH_GENERATORS__
#define __TESTH_GENERATORS__

#include "proc.h"

enum TestHGenerator { 
	TestH_EG		= -0x1,
	TestH_RF		= -0x2,
	TestH_RFM		= -0x3,
	TestH_RFT		= -0x4,
	TestH_RFTA		= -0x5,

	TestH_Gaussian	= 0x0,

	TestH_AR 		= 0x1,
	TestH_fBmSGA	= 0x2,
	TestH_3SPG		= 0x3,
	TestH_Hosk		= 0x4,
	TestH_Pax		= 0x5
};

typedef enum TestHGenerator gen;


void 	gen_PrintHeader	(gen g);
int 	gen_CheckGen	(gen g);

proc_Process* gen_GenProc (gen g, int n, double h);

proc_Process* gen_ExternGen (int N, char *name, double (*gen_func) (void), 
	tosig sign);


proc_Process* gen_ReadFile (const char *path, const char *name,
	tosig sig);
proc_Process* gen_ReadFileMatrix (const char *path, const char *name,
	tosig sig);
proc_Process* gen_ReadFileTime (const char *path, const char *name,
	tosig sig, int mult);

proc_Process* gen_Gaussian (int N, tosig sig);


proc_Process* gen_AggRenewal (int N, double h, int ren, int xm, tosig sig);


proc_Process* gen_fBmSequentialGenerationAlgorithm (int N, double h, 
	int scale_no, tosig sig);
proc_Process* gen_SimpleSelfSimilarProcessGenerator ();


proc_Process* gen_Hosking (int N, double h, tosig sig);
proc_Process* gen_Paxson  (int N, double h, tosig sig);

#endif