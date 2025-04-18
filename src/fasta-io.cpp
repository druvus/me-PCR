///////////////////////////////////////////////////////////////////
//
//		Multithreaded Electronic PCR (me-PCR) program
//
//		This implementation by:
//		Kevin Murphy
//		Children's Hospital of Philadelphia
//
//		Algorithm and most code by:
//		Gregory Schuler
//		Natonal Center for Biotechnology Information
//
//		Functions used by e-PCR to read FASTA files. 
//
///////////////////////////////////////////////////////////////////

// <2000-11-22 KPM: Many wonderful modifications
// 2002-11-22 KPM: Fixed destructors so all allocations are deallocated
// 2002-12-16 KPM: Terminated the parsed fasta buffer

#include <math.h>
#include <time.h>
#include "fasta-io.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

char chrValidAa[] = "ABCDEFGHIKLMNPQRSTUVWXZ-*";
char chrValidNt[] = "GATCNBDHKMRSVWY-";  


///////////////////////////////////////////////////////////////
//
//
//		class FastaSeq
//
//

FastaSeq::FastaSeq ()
{
	m_tag = NULL;
	m_def = NULL;
	m_seq = NULL; 
	m_len = 0;
}


FastaSeq::FastaSeq (char *def, char *seq)
{
	m_tag = NULL;
	m_def = NULL;
	m_seq = NULL; 
	SetDefline(def);
	SetSequence(seq);
	SetLength(strlen(seq));
}


FastaSeq::~FastaSeq ()
{
  // Deallocate m_tag, m_def, m_seq
  Clear(); 
}


const char * FastaSeq::Label() const
{
	static char zippo[] = "XXXXXX";
	return (m_tag==NULL) ? zippo : m_tag; 
}


const char * FastaSeq::Defline () const
{
	static char zippo[] = "XXXXXX  No definition line found";
	return (m_def==NULL) ? zippo : m_def; 
}


const char * FastaSeq::Title () const
{
	static char zippo[] = "No definition line found";
	if (m_def==NULL) return zippo;
	char *p;
	if ((p = strchr(m_def, ' ')))
	{
		while (*p && isspace(*p))  p++;
		if (*p) return p;
	}
	return m_def; 
}


const char * FastaSeq::Sequence () const
{
	return m_seq; 
}


int FastaSeq::Length() const
{
	return m_len;
}


// string is the first line of a sequence chunk, which contains "definition" information?
// m_def is a pointer to the definition line
// m_tag is an allocated string containing the identifier, e.g. L78833
void FastaSeq::SetDefline (char *string)
{
  if (m_tag) {
    MemDealloc(m_tag);
    m_tag = NULL;
  }
  if (m_def) {
    MemDealloc(m_def);
    m_def = NULL;
  }
  
  if (string) {
    
    int first_line_size = strcspn (string, "\r\n");
    if (!(MemAlloc(m_def, first_line_size+1))) {
      fprintf (stderr, "out of memory\n");
      exit (1);
    }
    strncpy (m_def, string, first_line_size);
    m_def[first_line_size] = '\0';
    
    char *p = m_def;
    // skip leading > and blanks ...
    while (*p && (*p == '>' || *p == ' ' || *p == '\t'))
      p++;
    if (*p) {
       int p_len = strlen(p);
       char *q = strpbrk(p, " \t");
       if (q)
	    p_len = q - p;
      if (!(MemAlloc(m_tag, p_len+1))) {
	fprintf (stderr, "out of memory\n");
	exit (1);
      }
       memcpy (m_tag, p, p_len);
       m_tag[p_len] = '\0';
    } else {
       // no tag found ...
       if (!(MemAlloc(m_tag,1))) {
	  fprintf (stderr, "out of memory\n");
	  exit (1);
       }
       m_tag[0] = '\0';
    }
  }
}


// m_seq points to the sequence data
void FastaSeq::SetSequence (char *string)
{
  if (string) {
    m_seq = string; 
  }
}

// m_len gets the length of the sequence data
void FastaSeq::SetLength (size_t len)
{
  m_len = len;
}


void FastaSeq::Clear()
{
  SetDefline(NULL); 
  SetSequence(NULL); 
  SetLength(0);
}




///////////////////////////////////////////////////////////////
//
//
//		class FastaFile
//
//


FastaFile::FastaFile ()
{
  FastaFile(0);
}


FastaFile::FastaFile (int seqtype)
{
  m_file = NULL;
  m_seqtype = seqtype;
  m_name = NULL;
  m_seqs = NULL;
  m_seq_base = NULL;
  m_numseqs = 0;
  m_maxseqs = 0;
}


FastaFile::~FastaFile ()
{
  if (m_file)  Close();
  // I guess I should be really using vector or something ...  Oh, well.
  if (m_seqs) {
    for (unsigned i=0; i<m_numseqs; i++) {
      delete m_seqs[i];
    }
    MemDealloc(m_seqs);
  }
  if (m_seq_base) {
    MemDealloc(m_seq_base);
  }
}

unsigned FastaFile::NumSeqs (void) 
{
  return m_numseqs;
}

FastaSeq **FastaFile::Seqs (void) 
{
  return m_seqs;
}

bool FastaFile::Open (const char *fname, const char *fmode)
{
	if (m_file != NULL)
	{
		PrintError("FastaFile::Open();  WARNING: file already open");
		return FALSE;
	}

	m_file = fopen(fname,fmode);
	if (m_file == NULL)
	{
		char buf[500];
		sprintf (buf, "FastaFile::Open();  ERROR: unable to open sequence file '%s'", fname);
		PrintError(buf);
		return FALSE;
	}

	m_name = fname;
	
	return TRUE;
}


bool FastaFile::Close ()
{
	if (m_file == NULL)
	{
		PrintError("FastaFile::Close();  WARNING: file already closed");
		return FALSE;
	}

	fclose(m_file);
	m_file = NULL;
	return TRUE;
}


#define CHUNK 100000

// 10/19/01 kpm: modified this routine to not free buf, since buf is now re-used to hold the 
// ... "parsed" (upcased and \n-filtered) sequence data.
// 4/8/04 kpm: now returns a vector of FastaSeq objects (FASTA can contain more than one sequence)
void FastaFile::Read (void)
{
  time_t start_time = 0;  // 0 just to avoid warning
  size_t FileSize;
  
  FileSize = ePCR_FileSize(m_name);
  
  if (!ePCR_quiet) {
    start_time = time(NULL);
    if (!FileSize) {
      fprintf (stderr, "Sequence file [%s] is empty!\n", 
	       m_name);
      goto Exit;
    }
    fprintf (stderr, "Reading sequence file ...\n");
  }
  
  if (!IsOpen())
    goto Exit;

  if (!MemAlloc(m_seq_base,FileSize+1)) {
    fprintf (stderr, "Out of memory in FastaFile::Read\n");
    goto Exit;
  }

  // read in one fell swoop
  if (fread (m_seq_base, FileSize, 1, m_file) != 1) {
    fprintf (stderr, "Error reading fasta file in FastaFile::Read\n");
    goto Exit;
  }

  m_seq_base[FileSize] = '\0';

  // sequences are parsed out and put in m_seqs property
  ParseText(m_seq_base, m_seqtype);
    
 Exit:
    
  if (!ePCR_quiet) {
    fprintf (stderr, "\t%3d %% done\n", 100);
    fprintf (stderr, "Elapsed time reading the sequence file: %lu seconds\n\n", (unsigned long) (time(NULL) - start_time));
  }

  return;
}


// Organize the sequences in the fasta file into a vector of FastaSeq
// objects (which involves stripping the control characters from the
// data, among other things).
void FastaFile::ParseText (char *text, const char *alphabet)
{
  time_t start_time = 0;  // 0 just to avoid warning
  unsigned char charmap[256];

#ifdef N_MODE
  unsigned char is_acgt[128];
#endif

  if (!ePCR_quiet) {
    start_time = time(NULL);
    fprintf (stderr, "Parsing the sequence file ...\n");
  }
  
  if (text==NULL || *text==0)
    return;

	
  // Build the map for filtering and upcasing the sequence data.
  memset (charmap, 0, sizeof(charmap));
  unsigned char c;
  while ((c=*alphabet++)) {
    charmap[toupper(c)] = toupper(c);
    charmap[tolower(c)] = toupper(c);
  }

#ifdef N_MODE
  // Build the map for testing ambiguity (post filter map)
  memset (is_acgt, 0, sizeof(is_acgt));
  is_acgt['A'] = is_acgt['C'] = is_acgt['G'] = is_acgt['T'] = 1;
#endif

  // Copy bases from p2 to p1, omitting newlines.  As new sequences
  // are encountered in the file (prefixed by '>'), create them and
  // add them to our vector of sequences.
  const char *p2 = text;   // source pointer
  char *p1;                // destination pointer
  char chr;
  
  // Loop over sequences (not characters)
  while (*p2) {

    if (*p2 != '>') {
      fprintf (stderr, "Error: expected '>'.  This version of e-PCR does not support more than one sequence in a file\n");
      exit (1);
    }
    
    if (m_numseqs == m_maxseqs) {
      // expand storage for sequences
      m_maxseqs += 100;
      if (!MemResize(m_seqs, m_maxseqs * sizeof(FastaSeq *))) {
	fprintf(stderr, "Out of memory trying to allocate room for %u sequences\n", m_maxseqs);
	exit(1);
      }
    }
    m_seqs[m_numseqs++] = new FastaSeq;

    size_t first_line_size = strcspn (p2, "\r\n");
    first_line_size += strspn (p2+first_line_size, "\r\n");

    // Record the description of this sequence (skipping the '>' character).
    m_seqs[m_numseqs-1]->SetDefline( (char *)p2 + 1 );   
    
    // Record the pointer to the start of the sequence and a fake sequence size (fix up later).
    const char *seq_start = p2+first_line_size;
    m_seqs[m_numseqs-1]->SetSequence ( (char *)seq_start );

    p2 += first_line_size;   // skip over description line
    p1 = (char *) p2;

#ifdef N_MODE
	int in_n_block = 0;
	int num_n_blocks_allocated = 0;
	int num_n_blocks = 0;
#endif

	// Loop over bases from the sequence.  Stop when we reach the
	// end of the buffer (zero-terminated) _or_ when we reach a
	// '>', signalling another sequence ('>' is not a valid
	// character in a FASTA file, outside of the description
	// line.)
	while ((chr = *p2) != '\0' && chr != '>') {
	  p2++;
	  if ((chr=charmap[chr])) {
#ifdef N_MODE
	    if (in_n_block) {
	      if (is_acgt[chr]) {
		in_n_block = 0;
	      }
	    } else {
	      if (!is_acgt[chr]) {
		if (num_n_blocks_allocated == num_n_blocks) {
		}
	      }
	    }
	      
#endif
	    *p1++ = chr;
	  }
	}  // end while
	
	// NOTA BENE: p2 is now pointing to either a '\0' or a '>'

	// sanity check
	if (chr == '>') {
	  if (*(p2-1) != '\n' && *(p2-1) != '\r') {
	    fprintf (stderr, "Error: unexpected '>' encountered not at the beginning of a line.\n");
	    exit(1);
	  }
	}
	
	// Terminate the newly moved buffer if need be.  This will
	// cause the search code to blow up if it tries to read past
	// the end of the sequence.  Note that we stomp on all the now
	// unused bytes to help turn up any problems reading past the
	// end of the sequence during the search phase.  If no
	// newlines were encountered during the parse, so that there
	// is no extra space following the (unmoved) data, the memset
	// does nothing (p2-p1 == 0).

	memset (p1, '\0', p2-p1);
	// harmless but unnecessary: *p1 = '\0';

	// Record the real length of the sequence
	m_seqs[m_numseqs-1]->SetLength( p1 - seq_start );

	if (!ePCR_quiet) {
		fprintf (stderr, "\t%3d %% done\n", 100);
		fprintf (stderr, "Elapsed time parsing the sequence file: %lu seconds\n\n", (unsigned long) (time(NULL) - start_time));
	}

  } // end while

  return;

} // end ParseText


void FastaFile::ParseText ( char *text, int seqtype)
{
  if (seqtype == SEQTYPE_AA)
    ParseText(text, chrValidAa);
  else if (seqtype == SEQTYPE_NT)
    ParseText(text,  chrValidNt);
  else {
    PrintError("ParseText(text,code);   ERROR: Invalid code");
    exit(1);
  }
  return;
} // end ParseText

