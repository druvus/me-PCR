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
//		Functions used by e-PCR to match STS primers against
//		sequences. 
//
///////////////////////////////////////////////////////////////////

// CHANGES:
// 11/19/02 Kevin Murphy: changed bzero() to memset() for portability.
// 11/20/02 KM: removed wordsize limit of 8
// 11/21/02 KM: now STS primers are mapped to uppercase;
//              fixed bug where large primers (>100) crashed the program.
// 11/22/02 KM: fixed nasty bug where the thread_args structure's num_hits
//              member wasn't getting initialized because of a typo!!!
// 11/22/02 KM: got rid of unused macro constants: WSIZE, ASIZE, etc.
// 12/05/02 KM: added ability to handle pcr size ranges, e.g. "200-220".
// 12/06/02 KM: now interprets ambiguous base characters correctly in STS's
//              -- but still no N's are allowed in the hash word!
// 12/09/02 KM: now hash values can occur at offsets other than
//              len-m_wsize in STS's.
// 12/10/02 KM: fixed reporting of right sequence position from (pos1+pos2) to (pos2+1).
// 12/12/02 KM: added handling for the I switch
// 12/16/02 KM: fixed IUPAC handling in the reverse direction
// 12/17/02 KM: fixed reads beyond the end of the sequence
// 12/19/02 KM: fixed uninitialized IUPAC mapping table
// 12/20/02 KM: fixed faulty discarding of redundant hits (when using multi threads)
// 01/10/02 KM: changed margin handling for STS's with a ranged size
// 07/28/03 KM: applied Peter Chines' patch and added '-' to mean "pcr size unknown".

#include <assert.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <ctype.h>

#ifndef __MWERKS__
#define _REENTRANT    /* do I need this? */
#include <pthread.h>
//#include <thread.h>

//kpm changed strings.h to string.h ok?
#include <string.h>
#endif

#include "stsmatch.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif


#define AMBIG 100

unsigned ePCR_STS_line_length = ePCR_MAX_STS_LINE_LENGTH_DEFAULT;

unsigned ePCR_default_pcr_size = ePCR_DEFAULT_PCR_SIZE_DEFAULT;

unsigned ePCR_iupac_mode = ePCR_IUPAC_MODE_DEFAULT;

char _scode[128];
int  _scode_inited;
char _compl[128];
int  _compl_inited;

// Our match algorithm for STS's containing ambiguous base characters
// will use the following lookup table.  The size (64K) seems extreme, but I am 
// going for speed over memory.  A reasonable alternative would be 
// to use 2^12 or 2^10 slots and do bit shifting and masking to create 
// the indices.  The table needs to be initialized by _init_IUPAC_match_matrix().
unsigned char *_IUPAC_match_matrix;
unsigned int _IUPAC_match_matrix_inited;
struct {
  char base;
  char *matches;
} _IUPAC_mapping[] = {
    {'A', "A"},
    {'C', "C"},
    {'G', "G"},
    {'T', "TU"},
    {'U', "TU"},
    {'R', "AGR"},
    {'Y', "CTUY"},
    {'M', "ACM"},
    {'K', "GTUK"},
    {'S', "CGS"},
    {'W', "ATUW"},
    {'B', "CGTUYKSB"},
    {'D', "AGTURKWD"},
    {'H', "ACTUYMWH"},
    {'V', "ACGRMSV"},
    {'N', "ACGTURYMKSWBDHVN"},
    {'\0', ""}
};



void init_scode ();
void init_compl ();
inline char *reverse (const char *from, int len, char *to);

// Poor EBCDIC users ...
#define MY_TOUPPER(c) ((c)&~32)


// Initialize the IUPAC mapping matrix used by seqmcmp_ambig()

void init_IUPAC_match_matrix (void) {
  int i = 0;

  if (ePCR_iupac_mode) {

    if (!MemAlloc (_IUPAC_match_matrix, 65536)) {
      fprintf (stderr, "out of memory\n");
      exit (1);
    }
    memset (_IUPAC_match_matrix, '\0', 65536);

    while (_IUPAC_mapping[i].base != '\0') {
      char *s = _IUPAC_mapping[i].matches;
      while (*s) {
	// STS value (.base) is in lower byte ...
	_IUPAC_match_matrix[(((unsigned short)(*s)) << 8) + (unsigned char)_IUPAC_mapping[i].base] = 1;
	s++;
      }
      i++;
    }
    _IUPAC_match_matrix_inited = 1;
  }
}

void free_IUPAC_match_matrix (void) {
  if (_IUPAC_match_matrix_inited && _IUPAC_match_matrix) {
    MemDealloc (_IUPAC_match_matrix);
  }
}



// Create an array for quick determination of ambiguity
unsigned char _ambig[256];
int ambig_inited;
#define PRIMER1 1
#define PRIMER2 2
void init_ambig (void) {
  _ambig['B'] = 1;	/* B = C, G or T */
  _ambig['D'] = 1;	/* D = A, G or T */
  _ambig['H'] = 1;	/* H = A, C or T */
  _ambig['K'] = 1;	/* K = G or T */
  _ambig['M'] = 1;	/* M = A or C */
  _ambig['N'] = 1;	/* N = A, C, G or T */
  _ambig['R'] = 1;	/* R = A or G (purines) */
  _ambig['S'] = 1;	/* S = C or G */
  _ambig['V'] = 1;	/* V = A, C or G */
  _ambig['W'] = 1;	/* W = A or T */
  _ambig['X'] = 1;	/* X = A, C, G or T */
  _ambig['Y'] = 1;	/* Y = C or T (pyrimidines) */
}


inline char * new_String (const char *str)
{
	if (str && *str)
	{
		char *str2 = new char[1+strlen(str)];
		strcpy(str2,str);
		return str2;
	}
	return 0;
}

inline void delete_String (char *str)
{
	delete [] str;
}


int InitStsList (const char *fname);


///////////////////////////////////////////////////////////////
//
//
//		class PCRmachine
//
//



PCRmachine::PCRmachine ()
{
	init_scode();
	init_compl();
	init_IUPAC_match_matrix();
	init_ambig();
	m_file = NULL;
	SetWordSize(ePCR_WDSIZE_DEFAULT);
	SetMargin(ePCR_MARGIN_DEFAULT);
	SetMismatch(ePCR_MMATCH_DEFAULT);
	SetThreePrimeMatch(ePCR_THREE_PRIME_MATCH_DEFAULT);
	max_pcr_size = 0;
	m_sts_table = NULL;
	m_last_global_sts = NULL;
}


PCRmachine::~PCRmachine ()
{
  STS *s = m_last_global_sts;
  free_IUPAC_match_matrix();
  if (!ePCR_quiet) fprintf (stderr, "Deleting STS's\n");
  while (s) {
    STS *next_s = s->global_prev;
    delete s;
    s = next_s;
  }
  delete [] m_sts_table;
  if (m_file) {
    if (!ePCR_quiet) fprintf (stderr, "Closing STS file\n");
    fclose(m_file);
  }
}


void PCRmachine::SetWordSize (int wdsize)
{
  if (wdsize < ePCR_WDSIZE_MIN || wdsize > ePCR_WDSIZE_MAX) {
    fprintf (stderr, "Error: wordsize (W=) of %d must be between %d and %d, inclusive", 
	     (int) wdsize,
	     (int) ePCR_WDSIZE_MIN,
	     (int) ePCR_WDSIZE_MAX);
    exit(1);
  } else

    m_wsize = wdsize;

  // Allocate storage for the hash table

  unsigned int a, i;
  // Each added wordsize bit actually implies 2 added bits of data
  // (the two bits needed to encode A,C,G or T), so each added wordsize
  // bit multiplies the number of possible hash values by 4
  for (i=0, a=1; i<m_wsize; i++)  { a *= 4; }
  
  m_asize = a;
  m_mask = a-1;

  if (!ePCR_quiet) fprintf (stderr, "m_asize=%u, m_mask=0x%x\n", m_asize, m_mask);

}


int PCRmachine::GetWordSize (void)
{
	return m_wsize;
}


void PCRmachine::SetThreePrimeMatch (unsigned three_prime_match)
{
  if (three_prime_match < ePCR_THREE_PRIME_MATCH_MIN) {
    fprintf (stderr, "Error: three_prime_match (X=) of %u must be greater than or equal to %u",
	     three_prime_match,
	     (unsigned) ePCR_THREE_PRIME_MATCH_MIN);
    exit(1);
  } else

    m_three_prime_match = three_prime_match;

  if (!ePCR_quiet) fprintf (stderr, "m_three_prime_match=%u\n", m_three_prime_match);

}


unsigned PCRmachine::GetThreePrimeMatch (void)
{
	return m_three_prime_match;
}


void PCRmachine::SetMargin (int margin)
{
  if (margin < ePCR_MARGIN_MIN || margin > ePCR_MARGIN_MAX) {
    fprintf (stderr, "Error: margin (M=) of %d must be between %d and %d, inclusive", 
	     (int) margin,
	     (int) ePCR_MARGIN_MIN,
	     (int) ePCR_MARGIN_MAX);
  } else {
    m_margin = margin;
  }
}


int PCRmachine::GetMargin (void)
{
	return m_margin;
}


void PCRmachine::SetMismatch (int mmatch)
{
  if (mmatch < ePCR_MMATCH_MIN || mmatch > ePCR_MMATCH_MAX) {
    fprintf (stderr, "Error: mismatch (N=) of %d must be between %d and %d, inclusive", 
	     (int) mmatch,
	     (int) ePCR_MMATCH_MIN,
	     (int) ePCR_MMATCH_MAX);
  } else {
    m_mmatch = mmatch;
  }

}


int PCRmachine::GetMismatch (void)
{
	return m_mmatch;
}




int PCRmachine::ReportHit (char *line, const char *seq_label, int pos1, int pos2, const STS *sts, long offset)
{
	fseek(m_file, offset, SEEK_SET);
	if (fgets(line, ePCR_STS_line_length+2, m_file))
	{
	        // Chomp off the end-of-line
    	        char *eol = line + strlen(line) - 1;
		*(eol+1) = '\0';
		while (eol >= line && *eol == '\n' || *eol == '\r')
		  *eol-- = '\0';

		// Make 'line' point to a string containing just the STS ID
		// and 'p' point to the tail of the line following the PCR size
		// (the tail may be nonexistent).
		char *p = strchr(line,'\t');
		*p++ = 0;
		p = strchr(p+1,'\t');
		p = strchr(p+1,'\t');
		p = strchr(p+1,'\t');   // p points to a tab or is NULL
		
		ePCR_printf( "%s\t%d..%d\t%s%s\t(%c)\n",seq_label,pos1+1,pos2+1,line,(p?p:""),sts->direct);
	}
	else
	{
		fprintf(stderr, "Error reading file (from offset %ld)\n",offset);
	}
	return TRUE;
}






/*
 * We have 1-ePCR_threads worth of results.  Each thread's results are
 * in order of the offset at which the hit occurred.  There may be
 * overlap between adjacent threads' result sets.  We eliminate
 * overlap by throwing out any hits for a thread which occur within
 * the first overlap bytes, UNLESS the thread is the first one (has a
 * base offset of 0).  
 * 2002-11-20 KPM: added support for EPCR_STATS and also non-quiet 
 * reporting of total hits, which is needed for the threaded version.
*/
void PCRmachine::ReportHits (const char *seq_label, epcr_thread_args_t *a, int num_threads)
{
  int i;
  unsigned long hit;
  char *line;
  unsigned long hits = 0;
#ifdef EPCR_STATS
  unsigned long hash_hits = 0, hash_looks = 0, string_comparisons = 0;
#endif

  if (!MemAlloc (line, ePCR_STS_line_length+2)) {
    fprintf (stderr, "out of memory in ReportHits() (pretty tragic ;-)\n");
    exit (1);
  }

  for (i=0; i<num_threads; i++) {
#ifdef EPCR_STATS
    hash_hits += a[i].hash_hits;
    hash_looks += a[i].comparisons;
    string_comparisons += a[i].string_comparisons;
#endif
    for (hit=0; hit<a[i].num_hits; hit++) {
      if (a[i].offset != 0 && (size_t) a[i].hits[hit].pos2 < m_overlap) {
	// This _should be_ a redundant hit!
	if (!ePCR_quiet) {
	  fprintf (stderr, "skipping redundant hit at offset %lu\n", (unsigned long) (a[i].offset + a[i].hits[hit].pos1));
	}
      } else {
	// Report a non-redundant hit
	ReportHit (line, seq_label, a[i].offset + a[i].hits[hit].pos1, a[i].offset + a[i].hits[hit].pos2, a[i].hits[hit].sts, a[i].hits[hit].sts->m_offset);
	hits++;
      }
    }
  }

  MemDealloc (line);
#ifdef EPCR_STATS
  fprintf (stderr, "hash looks = %lu, hash hits = %lu, string comparisons = %lu\n", hash_looks, hash_hits, string_comparisons); 
  if (ePCR_quiet) {
    // if STATS is on, the user probably wants this number!
    fprintf (stderr, "Total hits = %lu\n", hits);
  }
#endif
  if (!ePCR_quiet) {
    fprintf (stderr, "Total hits = %lu\n", hits);
  }
}


void PCRmachine::RecordHit (epcr_thread_args_t *a, int pos1, int pos2, const STS *sts)
{
  if (!ePCR_quiet) {
    static long hits;
    hits++;
    fprintf (stderr, "Hit %ld, thread %d\n", hits, a->id);
  }

  if (a->num_hits == a->num_hits_allocated) {
    // allocate more potential hits
    a->num_hits_allocated += 20;
    if (!MemResize(a->hits, a->num_hits_allocated * sizeof(epcr_hit_t))) {
      fprintf (stderr, "Out of memory in PCRmachine::RecordHit\n");
      exit(1);
    }
  }
  a->hits[a->num_hits].pos1 = pos1;
  a->hits[a->num_hits].pos2 = pos2;
  a->hits[a->num_hits].sts = (STS *)sts;
  a->num_hits++;
}


void *PCRmachine::ThreadProc (void *args)
{
  // Process sequence database (FASTA format)
  PCRmachine *object_ptr = static_cast<PCRmachine *>(((epcr_thread_args_t *)args)->object_ptr);

  if (!ePCR_quiet) {
    size_t offset = ((epcr_thread_args_t *)args)->offset;
    size_t length = ((epcr_thread_args_t *)args)->length;
    int id = ((epcr_thread_args_t *)args)->id;
    fprintf (stderr, "Thread %d is about to start working at offset %d over %d bytes\n", id, (int)offset, (int)length);
  }
  
  object_ptr->ProcessSeqThread((epcr_thread_args_t *)args);

  return NULL;
}



// This routine divides up the task for the desired number of threads and
// starts them rolling.  

int PCRmachine::ProcessSeq (const char *seq_label, const char *seq_data, size_t seq_len)
{
  int i;
  epcr_thread_args_t *arg_array;
  pthread_t *threads;
  int num_threads = ePCR_threads;  // the canonical case; may be reduced, though

  // This is the minimum required data m_overlap between each
  // thread chunk.
  m_overlap = max_pcr_size + m_margin - 1;

  if (!ePCR_quiet)
    fprintf (stderr, "Processing seq: '%s': m_overlap is %lu (from max_pcr_size of %lu and m_margin of %lu)\n", 
	     seq_label,
	     (unsigned long)m_overlap, 
	     (unsigned long)max_pcr_size, 
	     (unsigned long) m_margin);
  
  if ((num_threads > 1) && ((num_threads+1)*m_overlap > seq_len)) {

    // This combination of threads, margin, max pcr size, and sequence
    // length is not valid, because each thread wants a chance to find
    // the largest possible STS in its allocated subchunk of the
    // sequence.  So let's step down the number of threads until we
    // find the minimum number that work.

    do {
      --num_threads;
    } while ((num_threads > 1) && ((num_threads+1)*m_overlap > seq_len));

    /*
    if (max_pcr_size > 10000)
      fprintf (stderr, "\tThe max pcr size looks suspiciously large.  You may want to edit your STS file.\n");
    if (ePCR_threads > 20)
      fprintf (stderr, "\tThe number of threads is quite large.  "
	       "It is not useful to use many more threads than you have processors.  Try fewer threads.\n");
    if (m_margin > 10000)
      fprintf (stderr, "\tThe margin looks suspiciously large.  Try a lower value if possible.\n");
    */

  }

  size_t chunk_size = (size_t) ceil((seq_len - (num_threads+1)*m_overlap)/(double)num_threads) + 2*m_overlap;

  if (m_overlap >= chunk_size && num_threads > 1) {
    fprintf (stderr, "Error: m_overlap of %lu exceeds chunk size of %lu.  ",
	     (unsigned long)m_overlap, (unsigned long)chunk_size);
    fprintf (stderr, "Try fewer threads.  Also, "); 
    fprintf (stderr, 
	     "This could be cause by an excessively large STS (max pcr size from file is %d, or a large margin (currently %d)\n",
	     (int) max_pcr_size, (int) m_margin);
    exit(1);
  }

  size_t offset = 0;
  size_t length = chunk_size;

  if (!MemAlloc (arg_array, num_threads*sizeof(epcr_thread_args_t))) {
    fprintf (stderr, "out of memory\n");
    exit (1);
  }
  memset (arg_array, '\0', sizeof(arg_array));

  size_t last_offset = offset-1;
  
  if (!MemAlloc (threads, num_threads*sizeof(pthread_t))) {
    fprintf (stderr, "out of memory\n");
    exit (1);
  }
  memset (threads, '\0', sizeof(threads));

  if (!ePCR_quiet)
    fprintf (stderr, "sequence length=%lu\n", (unsigned long) seq_len);

  for (i = 0; i<num_threads; i++) {
    arg_array[i].id = i;
    arg_array[i].object_ptr = this;
    arg_array[i].offset = offset;
    arg_array[i].data = (char *)seq_data + offset;
    arg_array[i].length = (i < num_threads - 1) ? length : seq_len - offset;
    arg_array[i].num_hits = 0;
    arg_array[i].num_hits_allocated = 0;
    arg_array[i].hits = NULL;
    int rv;

    if (!ePCR_quiet) {
      fprintf (stderr, "thread %d will search from offset = %lu to %lu (length = %lu)\n", 
	       (int)i,
	       (unsigned long)arg_array[i].offset,
	       (unsigned long)arg_array[i].offset + (unsigned long)arg_array[i].length - 1,
	       (unsigned long)arg_array[i].length);
      fprintf (stderr, "--> leading m_overlap = %u\n", (unsigned) (last_offset - arg_array[i].offset + 1)); 
      last_offset = (long)arg_array[i].offset + (long)arg_array[i].length - 1;
    }
    offset += (length - m_overlap);
    
    //    rv = pthread_create(threads+i, NULL, ((void *)(*)(void *))ThreadProc, (void *)(arg_array+i));

    rv = pthread_create(threads+i, NULL, ThreadProc, arg_array+i);

    if (rv != 0) {
      fprintf (stderr, "error starting thread %d of %d. Code=%d\n", i+1, num_threads,rv);
      exit(1);
    }

  } // end for

  // Now just wait for the threads to complete ...

  for (i = 0; i<num_threads; i++) {
    /* joins the threads */
    pthread_join(threads[i], NULL);
  }
  
  // Report on the results.  All the threads have deposited their results in their 
  // private args structures.  There may be some duplicate hits.
  ReportHits (seq_label, arg_array, num_threads);

  for (i = 0; i<num_threads; i++) {
    /* joins the threads */
    if (arg_array[i].hits) 
      MemDealloc (arg_array[i].hits);
  }
  
  if (!ePCR_quiet) fprintf (stderr, "after joining threads\n");

  MemDealloc(threads);
  MemDealloc(arg_array);
     
  return 0;  

} // end ProcessSeq



/* This is the actual search algorithm
 * seq_data is upcased, whitespace-stripped sequence data
 * _scode is an array of 128 bytes, with ACGT mapped to 0,1,2,3, and everything else
 * mapped to 100 (AMBIG).
 */
int PCRmachine::ProcessSeqThread (epcr_thread_args_t *args)
{
  size_t seq_len = args->length;
  const char * seq_data = args->data;
  int count = 0;
  time_t start_time = 0;  // 0 just to suppress warning ...
#ifdef EPCR_STATS
  args->hash_hits = 0;
  args->comparisons = 0;
  args->string_comparisons = 0;
#endif

#ifdef TIME_TRIAL
  clock_t start_clock;
  start_clock = clock();
#endif
  
  if (!ePCR_quiet) {
    start_time = time(NULL);
    fprintf (stderr, "Processing the sequence ...\n");
  }
  
  if (seq_data && seq_len > m_wsize)
    {
      unsigned int h;
      const char *p = seq_data;
      int i, j, k, pos, N;
      
      /* kpm: A nice simple hash.  Actually, we're just
       * discarding the unused bits from each character and squeezing 
       * as many 2-bit symbols as possible into a variable.
       */
      for(i=N=0, h=0; (unsigned)i<m_wsize; i++)
	{
	  /// Initialize the hash value h with the first 
	  h <<= 2;
	  if ((j=_scode[*p++]) ==AMBIG)
	    {
	      N = m_wsize;
	    }
	  else
	    {
	      if (N >0) N--;
	      h |= (unsigned int) j;
	    }
	}
      
      // Notice that pos is not involved in loop termination
      // and that throughout the loop, 
      // pos = (p-m_wsize) - seq_data
      // i.e., pos is "m_wsize behind p".
      // First time through, pos=0 and p=seq_data+m_wsize.

      for (pos=0; (size_t)(p-seq_data)<seq_len; ++pos)
	{
	  // If N > 0, it means there was an N within the last m_wsize
	  // characters, so we know we don't have a valid hash value
	  // to test.
	  if (N == 0)
	    {
	      STS *sts = m_sts_table[h];
	      while (sts)
		{ 
#ifdef DEBUG
		  fprintf (stderr, "hash hit: %s/%s\n", sts->pcr_p1, sts->pcr_p2);
#endif
#ifdef EPCR_STATS
		  args->hash_hits++;
#endif
		  /* 
		   * When a hash match occurs, p points to the
		   * character after the last character of the
		   * m_wsize-sized match region, and pos indexes the
		   * first character of the match region.
		   *
		   * k indexes the first character of the primer.
		   *
		   * Note: the comparison against 0 is rather
		   * regrettable since it only applies to the first
		   * handful out of millions of comparisons.  For real
		   * speed we should prime the algorithm up to the
		   * point where this comparison is not needed.
		   */
		  k = pos - sts->hash_offset;
		  if (k>=0)
		     count += Match(
				    seq_data+k,
				    seq_len-k,
				    k,
				    sts, 
				    args
				    );
		  sts = sts->next;
		}  // end while
#ifdef EPCR_STATS
	      args->comparisons++;
#endif
	    }  // end if N==0

	  // Update the hash value.  If an ambiguous base (e.g. "N")
	  // is encountered, the hash is disqualified for wordsize
	  // bases after that.  We don't automatically forward the
	  // pointers by wordsize, but it doesn't seem to affect the
	  // speed.
	  h <<= 2;
	  h &= m_mask;
	  if ((j=_scode[*p++]) == AMBIG)
	    { 
	      // Interestingly, it doesn't really pay to try to scan for an entire block of N's here.
	      // And it will only pay less and less as time goes on ... so we won't.
	      N = m_wsize;
	    }
	  else
	    {
	       if (N>0) N--;
	       h |= (unsigned int) j;
	    }
	  
#ifdef TIME_TRIAL	  
			// For the Mac, allow e-PCR to operate in the background and play nice 
#ifdef __MWERKS__
	  TimeSlice();
#endif			
	  
	  // Output a progress indication
	  if (!ePCR_quiet) {
	    static int percent_done;
	    static clock_t last_clock;
	    clock_t this_clock;
	    float fract;
	    if ((this_clock = clock()) != last_clock) {
	      last_clock = this_clock;
	      fract = (float)(p-seq_data)/seq_len;
	      if (fract > percent_done/100.0) {
		percent_done  = (int) ceil(fract*100);
		fprintf (stderr, "\t%3d %% done; time spent: %d secs; hash hits: %lu; hash looks: %lu\n", 
			 percent_done-1, (int) (time(NULL)-start_time), args->hash_hits, args->comparisons);
	      }
	    }
	  }
#endif // TIME_TRIAL	  
	}
    }
  

#ifdef TIME_TRIAL
  fprintf (stderr, "Elapsed time processing the sequence: %f seconds\n\n", 
	   (float) (clock()-start_clock) / CLOCKS_PER_SEC);
  fprintf (stderr, "\thash hits: %lu; hash looks: %lu; total sequence length: %lu\n", 
			 args->hash_hits, args->comparisons, seq_len);
#endif
  
  if (!ePCR_quiet) {
    fprintf (stderr, "\t%3d %% done\n", 100);
    fprintf (stderr, "Elapsed time processing the sequence: %lu seconds\n\n", (unsigned long)(time(NULL) - start_time));
  }
  return count;
}


/*
  PCRmachine::Match()

  Once a hash hit on the left primer has occurred, Match() determines
  via brute force whether the whole STS matches the underlying
  sequence.  In doing so, Match() uses the expected PCR size and
  allowable margin (M) to find the right primer.  In the worst case,
  Match() will perform 2*M+1 string comparisons (via seqmcmp()) trying
  to find the right primer.

  Note: Match() always starts its search at the expected PCR size and
  gradually widens its search from there.  This makes sense since the
  graph of found versus expected size is a steep bell curve centered
  on the expected size.

  Note: There is various annoying logic here which prevents Match()
  from searching before the beginning of left primer (and hence,
  indirectly, before the beginning of the entire sequence, yikes) and
  also from searching beyond the end of the entire sequence.

  Hits are recorded by calling RecordHit() as they are encountered.

 Return value: the number of hits.

 */
inline int PCRmachine::Match (
			      const char *seq,  // Pointer to the beginning of the left primer in the sequence
			      size_t seq_len,   // Length of entire remaining sequence
			      int k,            // k indexes the first character of the primer
			      const STS *sts,   // The STS we're looking for
			      epcr_thread_args_t *args  // For storing statistics
			      ) 
{
   size_t len_p1 = sts->p1_len;
   size_t margin = sts->margin;
   int count = 0;
   
#ifdef EPCR_STATS
   args->string_comparisons++;
#endif

  if ( ((sts->ambig_primer&PRIMER1)?seqmcmp_ambig(seq,sts->pcr_p1,len_p1,+1):seqmcmp(seq,sts->pcr_p1,len_p1,+1))==0)
    {
      size_t len_p2 = sts->p2_len;
      size_t lo_margin, hi_margin;
      size_t exp_size = sts->pcr_size;

      // boundary check: we're not allowed to start looking after the end of the sequence!
      if (exp_size > seq_len) {

	// We'd like to coerce exp_size, but let's see if we can even do that ...

	if (seq_len < len_p1 + len_p2)
	  // no, this STS can't possibly fit this close to the end of the sequence
	  return 0;

	// yes, set exp_size to the maximum possible value
	exp_size = seq_len;
	hi_margin = 0;
	
      } else {

	// Keep exp_size the same, and ...
	hi_margin = margin;
	// Make sure that hi_margin will not extend our search beyond the sequence end.
	if ( (hi_margin + exp_size) > seq_len)
	  hi_margin = seq_len - exp_size;   // Note that in this if clause we are guaranteed that seq_len >= exp_size
      }

      // assert: seq_len >= exp_size

      lo_margin = margin;
      if (lo_margin > exp_size - len_p1 - len_p2)
	lo_margin = exp_size - len_p1 - len_p2;  
      
      // assert: lo_margin >= 0, because when the STS was created, exp_size was coerced to be >= len_p1 + len_p2

      assert ((int)seq_len >= (exp_size-lo_margin));

      
      /*      fprintf (stderr, "m,lo,hi=%d,%d,%d, p1=%s, lenp1,lenp2,stslen=%d,%d,%d, seq_len=%d\n",
	      margin, lo_margin, hi_margin, sts->pcr_p1, len_p1, len_p2, sts->pcr_size, (int)seq_len); */

#define P2_SEQMCMP(ptr,p2,len,strand) ((sts->ambig_primer&PRIMER2)?seqmcmp_ambig(ptr,p2,len,strand):seqmcmp(ptr,p2,len,strand))

#ifdef EPCR_STATS
  args->string_comparisons++;
#endif

      const char *p = seq + (exp_size - len_p2);
      if (seq_len>=exp_size && P2_SEQMCMP(p,sts->pcr_p2,len_p2,-1)==0)
	{
	   RecordHit(args, k, k+exp_size-1, sts);
	   count++;
	}
      
      size_t i;
      for (i=1; i<=margin; ++i)
	{
#ifdef EPCR_STATS
  args->string_comparisons++;
#endif
	  if (i<=lo_margin && P2_SEQMCMP(p-i,sts->pcr_p2,len_p2,-1)==0)
	     {
		RecordHit(args, k, k+exp_size-i-1, sts);
		count++;
	     }
#ifdef EPCR_STATS
	   args->string_comparisons++;
#endif
	  if (i<=hi_margin && P2_SEQMCMP(p+i,sts->pcr_p2,len_p2,-1)==0)
	     {
		RecordHit(args, k, k+exp_size+i-1, sts);
		count++;
	     }
	}
    }
   return count;
}


void PCRmachine::InsertSTS (STS *sts, unsigned hash)
{
  // Use hash value as index into array, insert item
#ifdef DEBUG
  fprintf(stderr,"Inserting STS: hash = %d (0x%04x), hash offset = %d, p1 = %s, p2 = %s, margin = %d, size = %d, ambig_primer = %d\n",
	  hash, hash, sts->hash_offset, sts->pcr_p1, sts->pcr_p2, sts->margin, sts->pcr_size, sts->ambig_primer);
#endif
  sts->next = m_sts_table[hash];
  m_sts_table[hash] = sts;
  m_sts_count++;
  if (m_last_global_sts)
    sts->global_prev = m_last_global_sts;
  m_last_global_sts = sts;
}




int PCRmachine::ReadStsFile (const char *fname)
{
  unsigned long FileSize = 0;  // 0 just to avoid warning ...
  time_t start_time = 0;  // 0 just to avoid warning ...

  max_pcr_size = 0;

  if (!ePCR_quiet) {
    FileSize = ePCR_FileSize(fname);
    if (FileSize == 0) {
      fprintf(stderr,"STS File [%s] empty!  (On Mac, make sure this isn't an alias)\n",fname);
      return 0;
    }
    start_time = time(NULL);
    fprintf (stderr, "Reading STS file ...\n");
  }

  m_file = fopen(fname,"rb");
  if (m_file == NULL)
    {
      fprintf(stderr,"Error: unable to open STS file: [%s]\n",fname);
      exit(1);
    }
  
  m_sts_table = new STS*[m_asize];
  memset((void*)m_sts_table,0,m_asize*sizeof(STS*));
  
  char *p;
  char *line;
  char *check_strtol;

  if (!MemAlloc (line, ePCR_STS_line_length+2)) {
    fprintf (stderr, "out of memory\n");
    exit (1);
  }

  STS *sts;
  char *pcr_p1, *pcr_p2;
  int len_p1, len_p2;
  long offset =0;
  int bad1=0, bad2=0, bad3=0;
  int line_no = 0;
  
  // Make sure we can detect a line buffer over-run on the first line read
  line[ePCR_STS_line_length] = '\0';
  
  // kpm: for each sts, 2 search targets are set up, with the first m_size characters hashed for quick searching
  while (fgets(line, ePCR_STS_line_length+2, m_file))
    {
      char *rev_p1 = NULL, *rev_p2 = NULL;
      int pcr_size =0;
      int margin_to_use = m_margin;  // This may change for STS's with a size range
      int hash_offset;
      char ambig_primer = 0, ambig_primer_rev = 0;  
      
      if (line[ePCR_STS_line_length] != '\0' && line[ePCR_STS_line_length] != '\n') {
	fprintf (stderr, 
		 "Error: the maximum STS file line length (not including line terminator(s)), %d, has been exceeded.\n", 
		 ePCR_STS_line_length);
	fprintf (stderr, 
		 "  Rerun e-PCR with S=<n> where <n> is the number output by the following command line:\n");
	fprintf (stderr,
		 "  perl -ne '$max=length($_)-1 if length($_) > $max; END{print \"$max\\n\";}' < %s\n", fname);
	exit(1);
      }
      line_no++;
      
      if (line[0] == '#') goto next_line;    // ignore comments
      if (line[0] == '\n') goto next_line;   // ignore blank lines

      if ((p = strchr(line,'\t')) ==NULL) {
	fprintf(stderr, "ERROR: bad STS file format (should be: idTABprimerTABprimerTABsize), file '%s', line # %d [%d]\n", fname, line_no, __LINE__);
	exit(1);
      }
      p++;
      pcr_p1 = p;
      if ((p = strchr(p,'\t')) ==NULL) {
	fprintf(stderr, "ERROR: bad STS file format (should be: idTABprimerTABprimerTABsize), file '%s', line # %d [%d]\n", fname, line_no, __LINE__);
	exit(1);
      }
      *p++ = 0;
      pcr_p2 = p;
      if ((p = strchr(p,'\t')) ==NULL) {
	fprintf(stderr, "ERROR: bad STS file format (should be: idTABprimerTABprimerTABsize), file '%s', line # %d [%d]\n", fname, line_no, __LINE__);
	exit(1);
      }
      *p++ = 0;

      if (*p == ' [%d]\n' || *p == '\0') {
	fprintf(stderr, "ERROR: bad STS file format (should be: idTABprimerTABprimerTABsize), file '%s', line # %d [%d]\n", fname, line_no, __LINE__);
	exit(1);
      }

      pcr_size = strtol(p, &check_strtol, 10);
      if (*check_strtol != '\n' && *check_strtol != '\0' && *check_strtol != '\t') {
	fprintf(stderr, "ERROR: size should be integer but is '%s'; bad STS file format (should be: idTABprimerTABprimerTABsize), file '%s', line # %d [%d]\n", p, fname, line_no, __LINE__);
	exit(1);
      }

      // Allow a pcr size range, in which case we use the average with a range (plus the normal margin).
      // This is a kludge - we should really use a lo/hi pair.
      if (pcr_size > 0) {
	while (*p && *p != '\t') {
	  if (*p == '-') {
	    int right = atoi(p+1);
	    int pcr_mod;
	    if (right == 0) {
	      fprintf (stderr, "Invalid PCR size value at line %d\n", line_no);
	      exit(1);
	    }
	    pcr_size += right;
	    pcr_mod = pcr_size % 2;
	    pcr_size /= 2;
	    margin_to_use += (right - pcr_size + 1);  // By this logic we will check 1 or 2 additional bytes
	    break;
	  }
	  p++;
	}
      } else {
	if (!( (*p == '-' && (!*(p+1) || isspace(*(p+1)))) || *p == '0')) {
	  fprintf (stderr, "Invalid PCR size value at line %d\n", line_no);
	  exit(1);
	}
      }

      if (pcr_size == 0) 
	pcr_size = ePCR_default_pcr_size;

      len_p1 = strlen(pcr_p1);
      len_p2 = strlen(pcr_p2);
      
      if (len_p1 + len_p2 > pcr_size) {
	// warn of an impossibly(?) small pcr_size, and coerce it
#ifdef DEBUG
	fprintf(stderr, "pcr size impossibly small at line %d of STS file: p1 = %s (len %d), p2 = %s (len %d), pcr_size = %d\n",
		__LINE__, pcr_p1, len_p1, pcr_p2, len_p2, pcr_size);
#endif
	bad3++;
	pcr_size = len_p1 + len_p2;
      }
      
      // kpm 2002-11-21: length of pcr_p2 wasn't being checked ...
      if (strlen(pcr_p1) < m_wsize || strlen(pcr_p2) < m_wsize)
	{
	  if (!ePCR_quiet)
	    fprintf (stderr, "\tWARNING [%s]: PCR primer shorter than word size \n",line);
	  bad1++;
	  offset = ftell(m_file);
	  goto next_line;
	}
      
      // 11/21/02 KPM: allowed arbitrarily large reverse primers
      if (!MemAlloc (rev_p1, len_p1+1)) {
	fprintf (stderr, "out of memory\n");
	exit (1);
      }
      if (!MemAlloc (rev_p2, len_p2+1)) {
	fprintf (stderr, "out of memory\n");
	exit (1);
      }
      
      for (char *p = pcr_p1; *p; p++) {
	*p = MY_TOUPPER(*p);
	if (ePCR_iupac_mode && _ambig[*p]) {
	  ambig_primer |= PRIMER1;
	  ambig_primer_rev |= PRIMER2;
	}
      }
      for (char *p = pcr_p2; *p; p++) {
	*p = MY_TOUPPER(*p);
	if (ePCR_iupac_mode && _ambig[*p]) {
	  ambig_primer |= PRIMER2;
	  ambig_primer_rev |= PRIMER1;
	}
      }

      reverse(pcr_p1,len_p1,rev_p1);		
      reverse(pcr_p2,len_p2,rev_p2);
      
      // Added for threads
      if (pcr_size > max_pcr_size) {
	if (!ePCR_quiet) {
	  fprintf (stderr, "old max=%d, new max=%d (LINENO=%d)\n", max_pcr_size, pcr_size, line_no);
	}
	max_pcr_size = pcr_size;
      }
      
      unsigned Hfor1, Hfor2;
      
      // Calculate the hash values for the primers.
      // Remember that the hash value is taken from the _end_
      // of the primer, not the beginning!
      // DO WE NEED TO CHECK ALL FOUR POSSIBILITIES???.
      // I'm getting rid of these two, since we don't even use their
      // hashes for anything.
      // !HashValue(rev_p2, Hrev2) ||
      // !HashValue(rev_p1, Hrev1) ||
      
      if ((hash_offset = HashValue(pcr_p1, len_p1, Hfor1)) == -1)
	{
	  if (!ePCR_quiet) fprintf(stderr,"can't make hash value for primer %s ...\n", line);
	  bad2++;
	  offset = ftell(m_file);
	  goto next_line;
	} 
      else
	{
	  sts = new STS(pcr_p1,rev_p2,'+',pcr_size,offset,margin_to_use,hash_offset,ambig_primer);
	  InsertSTS(sts,Hfor1);
	}
      
      if ((hash_offset = HashValue(pcr_p2, len_p2, Hfor2)) == -1)
	{
	  if (!ePCR_quiet) fprintf(stderr,"can't make hash value for %s ...\n", line);
	  bad2++;
	  offset = ftell(m_file);
	  goto next_line;
	}
      else
	{
	  sts = new STS(pcr_p2,rev_p1,'-',pcr_size,offset,margin_to_use,hash_offset,ambig_primer_rev);
	  InsertSTS(sts,Hfor2);
	}
      
      // future: think about putting STS's in an array for quick disposal?
      
      offset = ftell(m_file);
      
#ifdef __MWERKS__		
      // For a Mac, if priority is not max, give some cycles to the OS
      TimeSlice();
#endif
      
      // Output a progress indication
      if (!ePCR_quiet) {
	static int percent_done;
	static clock_t last_ticks;
	clock_t this_ticks;
	float fract;
	
	if ((this_ticks=clock()) > last_ticks) {
	  last_ticks = this_ticks;
	  fract = (float)offset/FileSize;
	  if (fract > percent_done/100.0) {
	    percent_done  = (int) ceil(fract*100);
	    fprintf (stderr, "\t%3d %% done\r", percent_done-1);
	  }
	}
      }
    next_line:
      
      // Make sure we can detect a too-long line on the next read
      line[ePCR_STS_line_length] = '\0';
      if (rev_p1) MemDealloc (rev_p1);
      if (rev_p2) MemDealloc (rev_p2);
    }
  
  
  /// NOTE: The file stays open!
  
  
  if (!ePCR_quiet) {
    fprintf (stderr, "\t%3d %% done\n", 100);
  }
  
  if (bad1)
    {
      fprintf(stderr,"\tWARNING: %d STSs have primer shorter than W (%d): not included in search ...\n", bad1, m_wsize);
    }
  if (bad2)
    {
      fprintf(stderr,"\tWARNING: %d primers have ambiguities which prevent computation of a hash value: not included in search ...\n", bad2);
    }
  if (bad3)
    {
      fprintf(stderr,"\tWARNING: %d STSs have a primer length sum greater than the pcr size: expected pcr size adjusted\n", bad3);
    }
  
  if (!ePCR_quiet) {
    fprintf (stderr, "Elapsed time reading the STS file: %lu seconds\n\n", (unsigned long) (time(NULL) - start_time));
  }
  
  MemDealloc (line);
  
  return 1;
}


// Compute a hash value for the specified primer.  Note that the hash
// value may not contain ambiguous bases (e.g. 'N').  If there is not
// a valid hash value at the end of the primer (i.e. the last m_wsize
// bases include an ambig), then the next earlier hash value is tried,
// until the beginning of the primer is reached.  If no valid hash
// value is found anywhere in the primer, -1 is returned.  Otherwise,
// the offset to the hash value is returned.

int PCRmachine::HashValue (const char *primer, int primer_len, unsigned &hash_value)
{
  unsigned int h;
  int i, j;
  const char *p;
  int offset = primer_len - m_wsize;
  
  do {    
    p = primer + offset;
    h = 0;
    // If i ever equals m_wsize, we have found a good hash value.
    for (i=0; (unsigned)i<m_wsize; ++i)
      {
	if ((j=_scode[*p++]) == AMBIG)
	  {
	    // Bad hash value means we have to jump back to the next
	    // possible hash value (we "add" one because of the offset
	    // decrement at the end of the while loop ....)
	    offset -= (m_wsize - i - 1);
	    break;
	  }
	h <<= 2;
	h |= (unsigned int) j;
      }  // endfor
    
    offset--;
  } while ( (offset >= 0) && ((unsigned)i<m_wsize) );
  
  if ((unsigned)i < m_wsize) {
    hash_value = 0x666;  // HEX value ;-)
    return -1;
  } else {
    hash_value = h;
    offset++; // bump it because we auto-decremented
    return offset;  
  }
} // endfunc


// Return 0 if two short pieces of sequence match, -1 otherwise,
// subject to m_mmatch (number of allowed mismatches) and m_three_prime_match
// (number of bases at the 3' end which MUST match).
// See also seqmcmp_ambig below. 
// Note that an 'N' in either primer or sequence can only match an 'N',
// so essentially, any STS with an N in it is only going to be found
// with N (mismatches) > 0.
// s2 is the STS string; s1 is the underlying sequence.
// strand is +1 if the 3' end is on the right and -1 if the 3' end is on the left.
int PCRmachine::seqmcmp (const char *s1, const char *s2, int len, int strand)
{
  const char *p1 = s1;
  const char *p2 = s2;
  int i, n;
  // Check mmatch so we can do slightly less work if the N parameter is 0
  if (m_mmatch != 0) {
    for (i=n=0; i<len; i++, p1++, p2++)
      {
	assert (*p1 != 0 && *p2 != 0);
	if (*p1 != *p2)
	  {
	    n++;
	    if (n>m_mmatch 
		|| (strand > 0 && i >= len-m_three_prime_match)
		|| (strand < 0 && i < m_three_prime_match)
		)
	      // No match if:
	      // *) we exceed the total allowed mismatches per primer, OR
	      // *) there is any mismatch with m_three_prime_match bases of 
	      // the 3' end of the primer.
	      return -1;
	  }
      }
    return 0;   // 0 means it matches (like strcmp)
  } else {
    for (i=0; i<len; i++, p1++, p2++)
      {
	assert (*p1 != 0 && *p2 != 0);
	if (*p1 != *p2)
	  return -1;
      }
    return 0;   // 0 means it matches (like strcmp)
  }
}



// Version of seqmcmp (above) that interprets ambiguous base symbols in the STS properly.
// s2 is the STS string; s1 is the underlying sequence.
int PCRmachine::seqmcmp_ambig (const char *s1, const char *s2, int len, int strand)
{
  const char *p1 = s1;
  const char *p2 = s2;
  int i, n;
  
  // Check mmatch so we can do slightly less work if the N parameter is 0
  if (m_mmatch != 0) {
    for (i=n=0; i<len; i++, p1++, p2++)
      {
	assert (*p1 != 0 && *p2 != 0);
	if (!_IUPAC_match_matrix[ ((unsigned short)(*p1) << 8) + *p2])
	  {
	    n++;
	    if (n>m_mmatch
		|| (strand > 0 && i >= len-m_three_prime_match)
		|| (strand < 0 && i < m_three_prime_match)
		)
	      // No match if:
	      // *) we exceed the total allowed mismatches per primer, OR
	      // *) there is any mismatch with m_three_prime_match bases of 
	      // the 3' end of the primer.
	      return -1;
	  }
      }
    return 0;   // 0 means it matches (like strcmp)
  } else {
    for (i=0; i<len; i++, p1++, p2++)
      {
	assert (*p1 != 0 && *p2 != 0);
	if (!_IUPAC_match_matrix[ ((unsigned short)(*p1) << 8) + *p2])
	  return -1;
      }
    return 0;   // 0 means it matches (like strcmp)
  }
}





///////////////////////////////////////////////////////////////
//
//
//		class STS
//
//



///// STS constructor and descructor

STS::STS (const char *p1, const char *p2, char d, int size, long offset, int margin_to_use, unsigned short p_hash_offset,char p_ambig_primer)
{
  if (*p1==0 || *p2==0)
    { 
      fprintf(stderr,"Error: empty primer: make sure the STS file is tab-delimited.  Byte offset in file is %ld\n",offset);
      exit(1);
    }
  
  next = NULL;
  
  direct = d;
  pcr_p1 = new_String(p1);
  p1_len = strlen(p1);
  pcr_p2 = new_String(p2);
  p2_len = strlen(p2);
  pcr_size = size;
  m_offset = offset;
  margin = margin_to_use;
  hash_offset = p_hash_offset;
  ambig_primer = p_ambig_primer;
  global_prev = NULL;
}

STS::~STS ()
{
  delete_String(pcr_p1);
  delete_String(pcr_p2);
}



/////////////////// Misc Utilities ///////////////////////

void init_scode ()
{
  if (!_scode_inited)
    {
      int i;
      for (i=0; (unsigned)i<sizeof _scode; ++i)
	_scode[i] = AMBIG;
      
      _scode['A'] = 0;
      _scode['C'] = 1;
      _scode['G'] = 2;
      _scode['T'] = 3;
      
      _scode_inited =1;
    }
}



void init_compl ()
{
  if (!_compl_inited)
    {
      _compl['A'] = 'T';
      _compl['C'] = 'G';
      _compl['G'] = 'C';
      _compl['T'] = 'A';
      
      /* ambiguity codes */
      _compl['B'] = 'V';	/* B = C, G or T */
      _compl['D'] = 'H';	/* D = A, G or T */
      _compl['H'] = 'D';	/* H = A, C or T */
      _compl['K'] = 'M';	/* K = G or T */
      _compl['M'] = 'K';	/* M = A or C */
      _compl['N'] = 'N';	/* N = A, C, G or T */
      _compl['R'] = 'Y';	/* R = A or G (purines) */
      _compl['S'] = 'S';	/* S = C or G */
      _compl['V'] = 'B';	/* V = A, C or G */
      _compl['W'] = 'W';	/* W = A or T */
      _compl['X'] = 'X';	/* X = A, C, G or T */
      _compl['Y'] = 'R';	/* Y = C or T (pyrimidines) */
      
      _compl_inited =1;
    }
}

inline char *reverse (const char *from, int len, char *to)
{
  const char *s;
  char *t;
  
  for (s = from+len-1, t = to; s >= from; --s, ++t)
    if ((*t = _compl[*s]) == 0)
      *t = 'N';
  *t = '\0';
  return to;
}


