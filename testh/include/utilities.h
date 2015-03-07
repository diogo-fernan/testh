#ifndef __TESTH_UTILITIES__
#define __TESTH_UTILITIES__

#include <time.h>


void* util_MemMalloc (size_t num);
void* util_MemCalloc (size_t num, size_t size);
void* util_MemRealloc (void *p, size_t num);
void* util_MemFree (void *p);

void util_MemWr (long int B);

void util_TimeIt (clock_t *c);
void util_TimeWr (clock_t *start);
void util_TimeReset (clock_t *c);

char* util_GetHostname ();
char* util_GetUsername ();
char* util_GetCwd ();

#endif