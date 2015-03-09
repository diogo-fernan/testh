#ifndef __TESTH_IO__
#define __TESTH_IO__

#include "gen.h"

#include <stdio.h>

#define ON 	1
#define OFF	0

#define OK	1
#define ERR -1

enum TestHVerbosity {
	TestH_NONE 		= 0x0, 
	TestH_LOW 		= 0x1, 
	TestH_MEDIUM 	= 0x2, 
	TestH_HIGH 		= 0x3
};

typedef enum TestHVerbosity verb;


extern verb TestHVerbosity;
extern int  TestHPrintPlain;
extern int  TestHPrintSep;
extern int  TestHPrintHeader;
extern int  TestHPrintMemCPU;
extern int  TestHEstPrintH;
extern int  TestHEstWrToFile;

// http://www.ibm.com/developerworks/linux/library/l-tip-prompt/
#define CGRAY_GRAY  "\e[30;1;47m"
#define CRED_ORANGE "\e[33;1;41m"
#define CGRAY_BLUE  "\e[34;1;47m"
#define CRESET      "\e[0m"

#define TESTH_printVar(var) \
	fprintf(stdout, "%s", #var);

	// Printouts
void io_PrintTestH ();
void io_PrintInit (const char *params, const char *argv);
void io_PrintDone ();
void io_PrintSep ();
void io_PrintError (const int err, const char *format, ...);
void io_PrintErr (const int err, const char *format, ...);

int io_CheckH (double h);
int io_CheckPowerOfTwo (int num);
int io_CheckPowerOf (int num, int exp);
// int io_PowerOfExp (int num, int exp);

	// Files
FILE* io_FileOpen (const char *path, const char *mode);
void io_FileClose (const char *path, FILE *f);
void io_FileClean (const char *path);
int io_FileLines (FILE *f);
int io_FileColumns (FILE *f);
int io_FileGetNum (FILE *f, double *p, gen g);
void io_FileWr (const char *path, const char *mode, 
	const double d1, const double d2);

#endif