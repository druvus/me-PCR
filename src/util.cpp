///////////////////////////////////////////////////////////////////
//
//		Electronic PCR (e-PCR) program
//
//              Original author:
//		Gregory Schuler
//		Natonal Center for Biotechnology Information
//
//              Extensively modified by:
//              Kevin Murphy
//              Children's Hospital of Philadelphia
//
//		Utility functions.
//
///////////////////////////////////////////////////////////////////
#define _REENTRANT    /* for threads */

#include <stdlib.h>
#include <errno.h>
#include <stdarg.h>

#ifdef __MWERKS__
#include <ctype.h>
#include <types.h>
#include <stat.h>
#include <events.h>
#include <SIOUX.h>
#else
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include "util.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif


unsigned long ePCR_hits;
unsigned ePCR_quiet = ePCR_QUIET_DEFAULT;
unsigned ePCR_priority = ePCR_PRIORITY_DEFAULT;
char *ePCR_outfile = ePCR_OUTFILE_DEFAULT;
int ePCR_threads = ePCR_THREADS_DEFAULT;

#if defined(__MWERKS__) | defined(__TURBOC__)

int strcasecmp (const char *str1, const char *str2)
{
  int c;
  while (*str1 && *str2) {
    c = toupper(*str1) - toupper(*str2);
    if (c)
      return c;
    str1++;
    str2++;
  }
  if (*str1)
    return +1;
  else if (*str2)
    return -1;
  else
    return 0;
}
#endif


#ifdef __MWERKS__
// If priority is not max, we will give up some time to the operating system and other applications.
// We give ourselves around 1/10th sec. of processing time for every (ePCR_PRIORITY_MAX - ePCR_priority)
// hundredths of a second we give the rest of the machine.  Thus, a priority of max priority of 30
// and a priority of 20 gives us the same as we give the rest of the machine.
void TimeSlice (void)
{
    EventRecord    theEvent;

#if 0
   static Int32 sLastTimeWeWaited = 0;// initalize to anything
   // or static unsigned long sLastTimeWeWaited = 0;
   
   if (TickCount() > sLastTimeWeWaited + 6)
   {
       WaitNextEvent(proper arguments );
      // handle the event, if necc
      // if (event.should be handled)
      //     MyHandleEvent();

      sLastTimeWeWaited = TickCount();
   }
#endif

	if (ePCR_priority < ePCR_PRIORITY_MAX) {
		static clock_t last_tenth;
		clock_t this_clock = clock();
		// If we accidentally put this is a tight loop, do a fast sanity check
		if (this_clock == last_tenth)
			return;
		// Only allow ourselves ca. 1/10th second of activity (assuming each line read takes less than 1/10th)
		if ( this_clock-last_tenth > .1*CLOCKS_PER_SEC ) {
			// Flush the event queue; very chivalrous
			while (GetNextEvent( everyEvent, &theEvent) )	
        		SIOUXHandleOneEvent(&theEvent);
            // Now WAIT for the next event, which will allow the system to be response to a user.
        	if (WaitNextEvent( everyEvent, &theEvent, ePCR_PRIORITY_MAX - ePCR_priority, 0L ) )
        		SIOUXHandleOneEvent(&theEvent);
            last_tenth = clock();
        }
	}
}
#endif


// Wrapper to abstract the main e-PCR output routine
int ePCR_printf(const char *fmt, ...)
{
         va_list args;
         char buf[1000];
        int i;
        FILE *f;
        static int virgin = 1;
        static char *open_mode = "a";
	static int using_stdout = 0;

	ePCR_hits++;

	if (virgin && strcasecmp(ePCR_outfile, "stdout")==0)
	  using_stdout = 1;

	if (using_stdout) 
	  f = stdout;
	else {
	  f = fopen (ePCR_outfile, open_mode);
	  if (!f) {
	    fprintf (stderr, "ERROR: can't open output file '%s'\n", ePCR_outfile);
	    exit(1);
	  }
	}
 		
	if (virgin) {
	  virgin = 0;
	  open_mode = "a";
	}

         va_start(args, fmt);
         i=vsprintf(buf,fmt,args);
         va_end(args);

         if (strlen(buf)>sizeof(buf)-1) {
         	fprintf (stderr, "ERROR: Output line exceed limit of %lu characters\n", sizeof(buf)-1);
         	exit(1);
         }

	 if (!ePCR_quiet)
	   fprintf (stderr, "\n\tHIT: %s", buf);

         if (fwrite (buf, strlen(buf), 1, f) != 1) {
	   fprintf (stderr, "ERROR: error writing to output file: %s\n", strerror(errno));
	   exit(1);
	 }

	 if (!using_stdout)
	   fclose (f);

         return i;
}


// Return the number of bytes in a file
unsigned long ePCR_FileSize (const char *fname) {
	struct stat a_stat;
		
	if (stat (fname, &a_stat) != 0) {
		fprintf (stderr, "Error: opening '%s': errno=%d\n",
				fname, errno);
		exit(1);
	} else {
		return a_stat.st_size;
	}	
}  // FileSize()


bool __MemAlloc (void **ptr, int size)
{
	if (size > 0)
	{
		*ptr = malloc(size);
		if (*ptr != NULL) {
		  return TRUE;
		}
		else {
		  *ptr = NULL;
		  return FALSE;
		}
	}
	fprintf( stderr, "error: __MemAlloc called with size == 0\n");
	exit(1);
}


bool __MemResize (void **ptr, int new_size)
{
	void * old_ptr = *ptr;

	if (old_ptr == NULL)
	  // fixed a bug here.  KPM 4/3/02
		return __MemAlloc(ptr,new_size);

	if (new_size == 0)
		return __MemDealloc(ptr);

	void *new_ptr = realloc(old_ptr,new_size);
	if (new_ptr == NULL)
	{
		return FALSE;
	}

	*ptr = new_ptr;
	return TRUE;
}


bool __MemDealloc (void **ptr)
{
	if (*ptr != NULL)
	{
		free(*ptr);
		*ptr = NULL;
		return TRUE;
	}
	return FALSE;
}


#ifdef GENERATE_HTML

void PrintError (const char *message)
{
	fprintf(stdout,"<P><B>%s</B>\n",message);
}


void FatalError (const char *message)
{
	PrintError(message);
	exit(0);
}

#else

void PrintError (const char *message)
{
	fprintf(stderr,"%s\n",message);
}


void FatalError (const char *message)
{
	PrintError(message);
	exit(1);
}

#endif
