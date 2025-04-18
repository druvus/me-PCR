#ifndef __util_h__
#define __util_h__

#define _REENTRANT    /* for threads */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#ifdef _AIX
#include <strings.h>
#endif

#ifndef bool
#define bool int
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE 
#define FALSE 0
#endif


//// Quiet flag (1=no extra messages; 0=verbose mode)
#define ePCR_QUIET_DEFAULT    1
#define ePCR_QUIET_MIN        0
#define ePCR_QUIET_MAX        1

//// Priority (how much ePCR will slow up the machine)
//// Note: this is not a linear scale.  29 and 30 are VERY different
#define ePCR_PRIORITY_DEFAULT 10
#define ePCR_PRIORITY_MIN     0
#define ePCR_PRIORITY_MAX     30

//// Default output filename
#define ePCR_OUTFILE_DEFAULT  "stdout"

//// Default number of threads
#define ePCR_THREADS_DEFAULT 1

extern unsigned ePCR_quiet;
extern unsigned ePCR_priority;
extern char *ePCR_outfile;
extern int ePCR_threads;

extern unsigned long ePCR_hits;  // count the number of hits

#if defined(__MWERKS__) | defined(__TURBOC__)
	int strcasecmp (const char *s1, const char *s2);
#endif

#ifdef __MWERKS__
	void TimeSlice (void);
#endif

bool __MemAlloc(void **ptr, int size);
bool __MemResize(void **ptr, int newsize);
bool __MemDealloc(void **ptr);

#define MemAlloc(x,y)     __MemAlloc((void**)&(x),(y))
#define MemResize(x,y)    __MemResize((void**)&(x),(y))
#define MemDealloc(x)     __MemDealloc((void**)&(x))

unsigned long ePCR_FileSize (const char *fname);

void PrintError(const char *message);
void FatalError(const char *message);

int ePCR_printf(const char*fmt, ...);

#endif
