
#include "util.h"
#include "io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include <time.h>
#include <ctype.h>
#include <sys/param.h>
#include <unistd.h>


// #ifdef WINDOWS
#if defined(__WIN32) || defined(__WIN64) || defined(__TOS_WIN__) || defined(__WINDOWS__)
	#include <direct.h>
	#define TESTH_GetCwd	_getcwd
	#define TESTH_Username	"USERNAME"
#else
// #elif defined(unix) || defined(__unix__) || defined(__unix)
	#include <unistd.h>
	#define TESTH_GetCwd	getcwd
	#define TESTH_Username	"USER"
#endif

#ifndef HOST_NAME_MAX
	#define HOST_NAME_MAX 255
#endif



void* util_MemMalloc (
	size_t num)
{
	if (num <= 0)
		io_PrintErr (ERR, "invalid size (<= 0) in"
			" util_MemMalloc");
	void *p;
	errno = 0;
	p = malloc (num);
	if (p == NULL)
		io_PrintErr (errno, "malloc failed in"
			" util_MemMalloc");
	return p;
}
void* util_MemCalloc (
	size_t num,
	size_t size)
{
	if (num <= 0 || size <= 0)
		io_PrintErr (ERR, "invalid num/size (<= 0) in"
			" util_MemCalloc");
	void *p;
	errno = 0;
	p = calloc (num, size);
	if (p == NULL)
		io_PrintErr (errno, "calloc failed in"
			" util_MemCalloc");
	return p;
}
void* util_MemRealloc (
	void 	*p,
	size_t 	num)
{
	if (num <= 0)
		io_PrintErr (ERR, "invalid size (<= 0) in"
				" util_MemRealloc");
	errno = 0;
	p = realloc (p, num);
	if (p == NULL)
		io_PrintErr (errno, "realloc failed in"
			" util_MemRealloc");
	return p;
}
void* util_MemFree (
	void *p)
{
	if (p != NULL) {
		free (p);
		p = NULL; // not working
	}
	return NULL;
}

static double util_MemMB (
	long int B)
{
	return (double) B / (1024 * 1024);
}
void util_MemWr (B)
	long int B;
{
	if (TestHVerbosity > TestH_NONE &&
		TestHPrintPlain == OFF && TestHPrintMemCPU == ON) {
		if (B < 0)
			io_PrintError (ERR, "negative bytes in"
				" util_MemWr");
		fprintf (stdout, "\n\t  Memory: %5.6lf MB\n",
			util_MemMB (B));
	}
}


void util_TimeIt (
	clock_t *c)
{
	errno = 0;
	if ((*c = clock ()) == (clock_t) -1 )
		io_PrintError (errno, "clock failed to determine clock time in"
			" util_TimeIt");
}
void util_TimeWr (
	clock_t *start)
{
	if (TestHVerbosity > TestH_NONE &&
		TestHPrintPlain == OFF && TestHPrintMemCPU == ON) {
		clock_t end;
		util_TimeIt (&end);
		if (*start != (clock_t) -1 && end != (clock_t) -1)
			fprintf (stdout, "\tCPU time: %5.6lf s\n\n",
				((double) (end - *start)) / CLOCKS_PER_SEC);
		else
			io_PrintError (ERR, "invalid clock_t structures (-1) in"
				" util_TimeWr");
		end = (clock_t) -1;
	}
	*start = (clock_t) -1;
}
void util_TimeReset (
	clock_t *c) {
	if (c != NULL) {
		*c = (clock_t) -1;
		util_TimeIt (c);
	}
}

void util_Copy (
	double 	*v,
	double	*v2,
	int 	size)
{	
	if (v == NULL || v2 != NULL || size <= 0)
		io_PrintError (ERR, "invalid parameters in"
			" util_Copy");
	int p;
	v2 = (double*) util_MemMalloc (size);
	for (p=0; p<size; p++) {
		v2[p] = v[p];
	}
}

void util_BubbleSort (
	double 	*v, 
	int 	size)
{
	if (v == NULL || size <= 0)
		io_PrintError (ERR, "invalid parameters in"
			" util_BubbleSort");

	int i, j;
	double p;
	for (i=0; i<size; i++) {
		for (j=1; j<size-i; j++) {
			if (v[j-1] > v[j]) {
				p = v[j-1];
				v[j-1] = v[j];
				v[j] = p;
			}
		}
	}
}


char* util_GetHostname () 
{
	errno = 0;
	char *host = (char*) malloc (HOST_NAME_MAX * sizeof (char));
	// host = NULL;
// #ifdef HAVE_UNISTD_H
// #ifdef G_OS_UNIX
	if (gethostname (host, HOST_NAME_MAX) == ERR) {
		io_PrintError (errno, "gethostname failed to determine hostname in"
			" util_GetHostname");
		host = NULL;
	}
// #else
// getenv 
/* #else
	io_PrintError (errno, "Unable to determine hostname in"
			" util_GetHostname");
#endif */
	return host;
}
char* util_GetUsername () 
{
	errno = 0;
	char *user = NULL;
#ifdef HAVE_UNISTD_H
	user = getlogin ();
#else
	user = getenv (TESTH_Username);
#endif
	if (user == NULL)
		io_PrintError (errno, "Unable to determine username in"
			" util_GetUsername");
	return user;
}
char* util_GetCwd () 
{
	errno = 0;
	char *cwd = (char*) malloc (MAXPATHLEN * sizeof (char));
	if (!TESTH_GetCwd (cwd, MAXPATHLEN)) {
		io_PrintError (errno, "getcwd failed to determine current working directory in"
			" util_GetCwd");
		return NULL;
	}
	cwd[MAXPATHLEN - 1] = '\0';
	return cwd;
}


