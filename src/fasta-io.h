#ifndef __fasta_io_h__
#define __fasta_io_h__

#include "util.h"

#define SEQLINE_LEN_DEFAULT 60   // default sequence characters per line
#define SEQLINE_LEN_MAX    120   // max sequence characters per line
#define SEQLINE_LEN_MIN     10   // min sequence characters per line

extern char chrValidAa[];
extern char chrValidNt[];

#define SEQTYPE_AA 1
#define SEQTYPE_NT 2


class FastaSeq
{
public:
	FastaSeq ();
	FastaSeq (char *def, char *seq);
	~FastaSeq ();

	const char * Label() const;
	const char * Title() const;
	const char * Defline() const;
	const char * Sequence() const;
	int Length() const;

	void Clear();

protected:
	char *m_tag;
	char *m_def;
	char *m_seq;
	int   m_len;
	void *m_bogus;

#ifdef N_MODE
	struct {
	  char *n_block_start;
	  char *n_block_end;
	} *n_blocks;
#endif
	void SetDefline (char *string);
	void SetSequence (char *string);
	void SetLength (size_t len);

	friend class FastaFile;
};


class FastaFile 
{
public:
	FastaFile();
	FastaFile(int seqtype);
	~FastaFile();

	bool Open (const char *fname, const char *fmode);
	bool Close ();

	bool IsOpen () const
		{ return (m_file ==NULL) ? 0 : 1; }

	void ParseText (char *text, int seqtype);
	void ParseText (char *text, const char *alphabet);

	void Read (void);
	bool Write (FastaSeq &seq);

	unsigned NumSeqs(void);
	FastaSeq **Seqs(void);

protected:
	FILE *m_file;
	char *m_seq_base;    // Pointer to sequence data buffer for fasta file
	int m_seqtype;
	const char *m_name;  // file name; just a pointer, note
	// The following 3 vars comprise a roll-your-own vector.  Should use STL later -kpm
	FastaSeq **m_seqs;   // ptr to array of ptrs to sequence objects
	unsigned m_numseqs;  // number of active sequences in array
	unsigned m_maxseqs;  // maximum number of slots allocated in array
};


extern int WriteSeqLines (FILE *fd, const char *seq, int len, int linelen=SEQLINE_LEN_DEFAULT);


#endif
