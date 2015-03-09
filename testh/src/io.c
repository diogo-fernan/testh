
#include "io.h"
#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <stdarg.h>
#include <errno.h>


#define TESTH_err (format, ...) \
	fprintf (stderr, "%s:%d: " format, __FILE__, __LINE__, __VA_ARGS__)

// printf ("%s", strerror (errno));
// cpp -dM /usr/include/errno.h | grep 'define E' | sort -n -k 3

verb TestHVerbosity	  = TestH_MEDIUM;
int  TestHPrintPlain  = OFF;
int  TestHPrintSep    = ON;
int  TestHPrintHeader = ON;
int  TestHPrintMemCPU = OFF;
int  TestHEstPrintH   = ON;
int  TestHEstWrToFile = OFF;

// another global parameter to turn on or off regression 
// at end of estimation, for example?
// more for other

	// Printouts
void io_PrintTestH ()
{
	if (TestHVerbosity > TestH_NONE && TestHPrintPlain == OFF) {
		io_PrintSep ();
		fprintf (stdout, 
			" TestH -- API Summary\n\n"
			" - proc.h (process)\n"
			"\t* struct proc_Process\n"
			"\t* struct proc_ScalesConfig\n"
			"\t* enum tosig\n"
				"\t\tTestH_fGn\n"
				"\t\tTestH_fBm\n"
			"\t* enum tos\n"
				"\t\tTestH_INC\n"
				"\t\tTestH_RAND\n"
				"\t\tTestH_POW\n"
			"\tproc_Create{Process/Scales[Time]/Points/ScalesConfig}\n"
			"\tproc_Delete{Process/Scales/Points}\n"
			"\tproc_Print{Process/Scales/Points}\n"
			"\tproc_Check{Process/Scales/Points/ScalesConfig}\n"
			"\tproc_PrintProcessStruct\n"
			"\tproc_SizeOfProcess\n"
			"\tproc_Normalize\n"
			"\tproc_FractionalGaussianNoise\n"
			"\tproc_FractionalBrownianMotion\n"
			" - gen.h (generators)\n"
			"\t* enum gen\n"
				"\t\tTestH_EG\n"
				"\t\tTestH_RF\n"
				"\t\tTestH_RFT\n"
				"\t\tTestH_Gauss\n"
				"\t\tTestH_AR\n"
				"\t\tTestH_fBmSGA\n"
				"\t\tTestH_3SPG\n"
				"\t\tTestH_Hosk\n"
			"\tgen_ExternGen\n"
			"\tgen_ReadFile[Time]\n"
			"\tgen_Gaussian\n"
			"\tgen_AggRenewal\n"
			"\tgen_fBmSequentialGenerationAlgorithm\n"
			"\tgen_SimpleSelfSimilarProcessGenerator\n"
			"\tgen_Hosking\n"
			" - est.h (estimators)\n"
			"\t* enum est\n"
				"\t\tTestH_RS\n"
				"\t\tTestH_VT\n"
				"\t\tTestH_AMT\n"
				"\t\tTestH_EBP\n"
			"\test_RescaledRangeStatistics\n"
			"\test_VarianceTime\n"
			"\test_AbsoluteMomentsTime\n"
			"\test_EmbeddedBranchingProcess\n"
			" - batt.h (batteries)\n"
			"\tbatt_Standard\n"
			"\n");
	}
}

void io_PrintInit (params, argv)
	const char *params;
	const char *argv;
{
	if (TestHVerbosity > TestH_NONE && TestHPrintPlain == OFF) {
		fprintf (stdout, "\n Usage: %s\n\n", params);
		fprintf (stdout, "\thost: %s\n", util_GetHostname ());
		fprintf (stdout, "\tuser: %s\n", util_GetUsername ());
		fprintf (stdout, "\t cwd: %s\n", util_GetCwd ());
		fprintf (stdout, "\t bin: %s\n\n", argv);
	}
}
void io_PrintDone () 
{
	if (TestHVerbosity > TestH_NONE && TestHPrintPlain == OFF)
		fprintf (stdout, "\n\n done!\n\n");
}
void io_PrintSep () 
{
	if (TestHVerbosity > TestH_NONE && 
		TestHPrintSep == ON && TestHPrintPlain == OFF)
		fprintf (stdout, "\n%s- - - - - - - - - - - - - - "
			"- - - - - - - - - - - - - -%s\n\n",
			CGRAY_GRAY, CRESET);
}
void io_PrintError (
	const int err,
	const char *format,
	...)
{
	// if (TestHVerbosity > TestH_NONE) {
		va_list args;
		va_start (args, format);
		fprintf (stderr, "\n[!] %sERROR%s %d: %s\n", 
			CRED_ORANGE, CRESET, err, strerror (err));
		fprintf (stderr, "%s:%d: ", __FILE__, __LINE__);
		vfprintf (stderr, format, args);
		va_end (args);
		fprintf (stderr, "\n\n");
	// }
}
void io_PrintErr (
	const int err,
	const char *format,
	...)
{
	// if (TestHVerbosity > TestH_NONE) {
		va_list args;
		va_start (args, format);
		fprintf (stderr, "\n[!] %sERROR%s %d: %s\n", 
			CRED_ORANGE, CRESET, err, strerror (err)); // perror
		fprintf (stderr, "%s:%d: ", __FILE__, __LINE__);
		vfprintf (stderr, format, args);
		va_end (args);
		fprintf (stderr, "\n\n");
	// }
	exit (EXIT_FAILURE);
}

void io_PrintBits (
	unsigned long x, 
	int k)
{
/* void bits (unsigned long long i) {
	fprintf (stdout, "\n\t%llu (%lu) = ", i, sizeof (i) * 8);
	int c = 0;
	// unsigned long long mask = 0x8000000000000000L; // 8 bytes
	unsigned long long mask = 0x80000000L; // 4 bytes
	for (; mask != 0x00;) {
    	fprintf (stdout, "%c", !!(mask & i) + '0');
    	mask >>= 1;
	    if (!(++c % 4))
	    	fprintf (stdout, " ");
	}
	fprintf (stdout, "\n\n");
} */

/*    int i, n = CHAR_BIT * sizeof (unsigned long);
   unsigned long mask = (unsigned long) 1 << (n - 1);
   int spaces;
   lebool flag = FALSE;

   if (k > 0) {
      spaces = k - n;
      for (i = 0; i < spaces; i++)
         printf (" ");
   }
   for (i = 0; i < n; i++) {
      if (x & mask) {
         printf ("1");
         flag = TRUE;
      } else if (flag)
         printf ("0");
      else
         printf (" ");
      mask >>= 1;
   }
   if (k < 0) {
      spaces = -k - n;
      for (i = 0; i < spaces; i++)
         printf (" ");
   } */
}

int io_CheckH (
	double h)
{
	return h < 0.0 || h >= 1.0 
		? ERR
		: OK; 
}

int io_CheckPowerOfTwo (num)
	int num;
{
	return (num != 0) && ((num & (num - 1)) == 0);        
}

int io_PowerOfExp (num, b)
	int num;
	int b;
{
	if (num <= 0 || b <= 0)
		io_PrintErr (ERR, "invalid parameters in"
			" io_PowerOfExp");
	return (int) (log (num) / log (b));
}

int io_CheckPowerOf (num, b)
	int num;
	int b;
{
	if (num <= 0 || b <= 0)
		io_PrintErr (ERR, "invalid parameters in"
			" io_CheckPowerOf");
	// return fabs (remainder (log (num), log (b))) == 0.0;
	double d = log (abs (num)) / log (abs (b));

	if ((num > 0 && b > 0) || (num < 0 && b < 0)) {
		if (d == (int) d)
			return OK;
		else
			return 0;
	} else if (num > 0 && b < 0) {
		if ((int) d % 2 == 0)
			return OK;
		else
			return 0;
	} else
		return 0;
}


FILE* io_FileOpen (path, mode)
	const char *path;
	const char *mode;
{
	FILE *f;
	errno = 0;
	if ((f = fopen (path, mode)) == NULL)
		io_PrintErr (errno, "fopen failed opening %s in"
			" io_OpenFile", path);
	return f;
}
void io_FileClose (path, f)
	const char *path;
	FILE *f;
{
	errno = 0;
	if (f != NULL && 
		fclose (f) == EOF)
			io_PrintErr (errno, "fclose failed closing %s in"
				" io_FileClose", path);
}
void io_FileClean (path) 
	const char *path;
{
	FILE *f = io_FileOpen (path, "w");
	io_FileClose (path, f);
}
int io_FileLines (f)
	FILE *f;
{
	if (f == NULL)
		io_PrintErr (ERR, "invalid FILE pointer (NULL) in"
			" io_FileLines");
	int l = 0, ch = 0;
	while (EOF != (ch = fgetc(f))) {
		if (ch == '\n') {
			l++;
		}
	}
	return l;
}

int io_FileGetNum (f, p)
	FILE *f;
	double *p;
{
	if (f == NULL || p == NULL)
		io_PrintErr (ERR, "invalid FILE or double pointer (NULL) in"
			" io_FileGetNum");

	int s, m, d, dd, ch;
	char *str = NULL;
	s = m = d = dd = ch = 0;
	
	while (EOF != (ch = fgetc(f))) {
		// fprintf (stdout, " %c %d\n", ch, s);
		if (ch == '-')
			m = 1;
		if ((ch == '.' || ch == ',') && d != 1) {
			if (s > 0)
				d = 1;
		}
		else {
			if (isdigit (ch)) {
				if (s == 0) {
					if (m == 1) {
						s = 2;
						str = (char*) realloc (str, s * sizeof (char));
						str[0] = '-';
						str[1] = ch;
					} else {
						str = (char*) realloc (str, ++s * sizeof (char));
						str[0] = ch;
					}
				} else if (d == 1 && dd == 0) {
					s += 2;
					str = (char*) realloc (str, s * sizeof (char));
					str[s-2] = '.';
					str[s-1] = ch;
					dd = 1;
				} else {
					str = (char*) realloc (str, ++s * sizeof (char));
					str[s-1] = ch;
				}
			} else {
				if (s == 0)
					continue;
				str = (char*) realloc (str, (s + 1) * sizeof (char));
				str[s] = '\0';
				*p = atof (str);
				free (str);
				return OK;
			}
		}
	}
	return EOF;
}
void io_FileWr (path, mode, d1, d2)
	const char *path;
	const char *mode;
	const double d1;
	const double d2;
{
	FILE *f = io_FileOpen (path, mode);	
	fprintf (f, "%lf %lf\n", d1, d2);
	io_FileClose (path, f);
}
