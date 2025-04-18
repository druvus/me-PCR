#ifndef __stsmatch_h__
#define __stsmatch_h__

#include "util.h"

// This can be increased using S=## on the command line
#define ePCR_MAX_STS_LINE_LENGTH_DEFAULT 1022   // So line buffer is 1024

extern unsigned ePCR_STS_line_length;
extern unsigned ePCR_default_pcr_size;
extern unsigned ePCR_iupac_mode;
extern unsigned ePCR_3prime_bases_must_match;

/*
Implications of e-PCR wordsize on the number of MB of RAM
consumed by the hash table:
Wdsize  Hash_entries    MB
1       4       0.0
2       16      0.0
3       64      0.0
4       256     0.0
5       1024    0.0
6       4096    0.0
7       16384   0.1
8       65536   0.2
9       262144  1.0
10      1048576 4.0
11      4194304 16.0
12      16777216        64.0
13      67108864        256.0
14      268435456       1024.0
15      1073741824      4096.0
16      4294967296      16384.0
This version of e-PCR supports up to 16 bits, which would 
take 16GB of RAM.  Hope that's not too limiting.
*/

//// Word size
#define ePCR_WDSIZE_DEFAULT   11
#define ePCR_WDSIZE_MIN       3
#define ePCR_WDSIZE_MAX       16

//// Number of mismatches allowed
#define ePCR_MMATCH_DEFAULT   0
#define ePCR_MMATCH_MIN       0
#define ePCR_MMATCH_MAX       10

//// Margin (allowed deviation in product size)
#define ePCR_MARGIN_DEFAULT   50
#define ePCR_MARGIN_MIN       0
#define ePCR_MARGIN_MAX       10000

////  Number of 3' bases which must match (both primers)
////  The word size parameter implicitly forces matches on the 3' side of the left primer.
#define ePCR_THREE_PRIME_MATCH_MIN 0
#define ePCR_THREE_PRIME_MATCH_DEFAULT 1
//#define ePCR_THREE_PRIME_MATCH_MAX 1



//// PCR size to use if the stated size is 0.
//// Note that if the user specifies an impossibly low one, 
//// the program will issue a warning and coerce the size for
//// each STS to the sum of the primer lengths.
#define ePCR_DEFAULT_PCR_SIZE_DEFAULT 240
#define ePCR_DEFAULT_PCR_SIZE_MIN 1
#define ePCR_DEFAULT_PCR_SIZE_MAX 10000


#define ePCR_IUPAC_MODE_DEFAULT 0
#define ePCR_IUPAC_MODE_MIN 0
#define ePCR_IUPAC_MODE_MAX 1


class STS
{
	friend class PCRmachine;

public:
	STS  *next;      // pointer to next element in the linked list
	char *pcr_p1;    // left primer
	unsigned short   p2_len;    // length of left primer
	char *pcr_p2;    // right primer
	unsigned short   p1_len;    // length of right primer
	int   pcr_size;  // size of PCR amplicon
	int   margin;    // Margin to use when searching (
	unsigned short hash_offset;  // offset of the hash from the normal position
	char  ambig_primer;  // PRIMER1, PRIMER2, PRIMER1|PRIMER2, or 0
	char  direct;    // 'p' for plus, 'm' for minus
	long  m_offset;  // offset into STS primer file (beginning of line)
	STS *global_prev;  // allows deallocating all STS's quickly

	STS (const char *p1, const char *p2, char d, int size, long offset, int p_pcr_size_margin, unsigned short p_hash_offset,char p_ambig_primer);
	~STS ();

	//int Match (const char *seq, int margin, int &size) const;
};



typedef struct {
  int pos1;
  int pos2;
  STS *sts;
} epcr_hit_t;


typedef struct {
  int id;
  void *object_ptr;
  char * data;
  size_t offset;
  size_t length;
  epcr_hit_t *hits;
  unsigned long num_hits;
  unsigned long num_hits_allocated;
#ifdef EPCR_STATS  
  unsigned long hash_hits;
  unsigned long comparisons;
  unsigned long string_comparisons;
#endif
} epcr_thread_args_t;


class PCRmachine
{
public:
        static const size_t MIN_FILESIZE_FOR_THREADING = 100000;

	PCRmachine();
	virtual ~PCRmachine();

	int ReadStsFile (const char *fname);
	int ProcessSeqThread (epcr_thread_args_t *args);
	int ProcessSeq (const char *seq_label, const char *seq_data, size_t seq_len);
	static void *ThreadProc (void *args);

	void SetWordSize (int wdsize);
	int GetWordSize (void);
	void SetMargin (int margin);
	int GetMargin (void);
	void SetMismatch (int mmatch);
	int GetMismatch (void);
	void SetThreePrimeMatch (unsigned bases);
	unsigned GetThreePrimeMatch (void);
	unsigned long SizeStsFile (const char *fname);

	int   max_pcr_size;

	// Override the ReportHit() function for customized output of 
	// results.  Return value: FALSE to abort search, TRUE to continue

	virtual int ReportHit (
	        char *line,                 // space to read in the original STS file line
		const char *seq_label,      // Label for sequence
		int pos1, int pos2,         // STS endpoints, zero-based
		const STS *sts,             // STS that was hit
		long offset );              // Offset into STS file

protected:
	FILE *m_file;	// STS primer file
	STS **m_sts_table;
	unsigned long  m_sts_count;
	int   m_margin;
	int   m_mmatch;
	size_t m_overlap;
	unsigned int m_wsize;
	unsigned int m_asize;
	unsigned int m_mask;
	unsigned int m_three_prime_match;
	STS *m_last_global_sts;   // Pointer to chain of all STS's for convenient destruction

	void InsertSTS (STS *sts, unsigned hash);
	int HashValue (const char *primer, int primer_len, unsigned &hash);
	inline int Match (
			  const char *seq,
			  size_t seq_len, 
			  int k,
			  const STS *sts,
			  epcr_thread_args_t *args
        );
	inline int seqmcmp (const char *s1, const char *s2, int len, int strand);
	inline int seqmcmp_ambig (const char *s1, const char *s2, int len, int strand);
	void ReportHits (const char *seq_label, epcr_thread_args_t *a, int num_threads);
	void RecordHit (epcr_thread_args_t *a, int pos1, int pos2, const STS *sts);
};


#endif

/* Compilers disagree where these class constants should be initialized */
/* const size_t PCRmachine::MIN_FILESIZE_FOR_THREADING = 100000;*/
