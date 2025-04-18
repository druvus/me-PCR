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
// $Id: me-PCR.cpp,v 1.12 2008/02/01 03:27:20 murphy Exp $
///////////////////////////////////////////////////////////////////


// To test, try me-PCR test/test.sts test/test.fa

#include "stsmatch.h"
#include "fasta-io.h"

#ifdef __MWERKS__
#include <string.h>
#include "console.h"
#endif

#ifdef DMALLOC
#include "dmalloc.h"
#endif

// Modify this when you release a new version.
// First modify this.  Then cvs commit
// Then cvs tag, then 'make distrib' from 
// parent directory.
const char *release_version = "1.0.6";

static int Usage(void)
{
	fprintf(stderr,"\nme-PCR: Multithreaded Electronic PCR\n");
	fprintf(stderr,"Original author: Gregory Schuler, NCBI\n");
	fprintf(stderr,"Modified by Kevin Murphy, Children's Hospital of Philadelphia (murphy@genome.chop.edu)\n");
	fprintf(stderr, "Release %s\n", release_version);

	fprintf(stderr,"USAGE:  me-PCR stsfile seqfile [options]\n");
	fprintf(stderr,"OPTIONS:\n");
	fprintf(stderr,"\tM=##     Margin (default %d)\n",ePCR_MARGIN_DEFAULT);
	fprintf(stderr,"\tN=##     Number of mismatches allowed (default %d)\n",ePCR_MMATCH_DEFAULT);
	fprintf(stderr,"\tX=##     Number of 3' bases which must match (both primers) (default %d)\n",
		ePCR_THREE_PRIME_MATCH_DEFAULT);
	fprintf(stderr,"\tW=##     Word size (default %d)\n",ePCR_WDSIZE_DEFAULT);
	fprintf(stderr,"\tT=##     Number of threads (default 1)\n");
	fprintf(stderr,"\tO=file   Output file name (default %s)\n",ePCR_OUTFILE_DEFAULT);
	fprintf(stderr,"\tQ=##     Quiet flag\n");
	fprintf(stderr,"\t            0 = verbose progress messages\n");
	fprintf(stderr,"\t            1 = no progress messages (default)\n");
	fprintf(stderr,"\tS=##     Max. line length for the STS file (not counting line terminators) (default %d)\n",
		ePCR_MAX_STS_LINE_LENGTH_DEFAULT);
	fprintf(stderr,"\tZ=##     Default PCR size (default %d)\n",ePCR_DEFAULT_PCR_SIZE_DEFAULT);
	fprintf(stderr,"\tI=#      IUPAC flag\n");
	fprintf(stderr,"\t            0 = do not honor IUPAC ambiguity symbols in STS's (default)\n");
	fprintf(stderr,"\t            1 = honor IUPAC ambiguity symbols in STS's\n");


#ifdef __MWERKS__
	fprintf(stderr,"\tP=##     Priority (for Mac OS <10): 0-30 (default %d)\n",ePCR_PRIORITY_DEFAULT);
	fprintf(stderr,"\t            0 = takes longest\n");
	fprintf(stderr,"\t           30 = fastest, but completely hogs machine\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"To break out of program, press Control-Option-Command-Esc\n");
#endif
	return 1;
}



int main (int argc, char **argv)
{
	//// Gather arguments

	const char *stsfile = NULL;
	const char *seqfile = NULL;
	int margin = ePCR_MARGIN_DEFAULT;
	int mmatch = ePCR_MMATCH_DEFAULT;
	int wdsize = ePCR_WDSIZE_DEFAULT;
	unsigned three_prime_match = ePCR_THREE_PRIME_MATCH_DEFAULT;

#ifdef __MWERKS__
    argc = ccommand(&argv);
#endif

	int i;
	for (i=1; i<argc; ++i)
	{
		if (argv[i][1] == '=')         // X=value
		{
			if (argv[i][2] == 0)
				fprintf(stderr,"Missing value for %s\n",argv[i]);
			else if (argv[i][0] == 'M')
				margin = atoi(argv[i]+2);
			else if (argv[i][0] == 'N')
				mmatch = atoi(argv[i]+2);
			else if (argv[i][0] == 'W')
				wdsize = atoi(argv[i]+2);
			else if (argv[i][0] == 'O')
				ePCR_outfile = argv[i]+2;
			else if (argv[i][0] == 'Q')
				ePCR_quiet = atoi(argv[i]+2);
			else if (argv[i][0] == 'P')
				ePCR_priority = atoi(argv[i]+2);
			else if (argv[i][0] == 'S')
				ePCR_STS_line_length = atoi(argv[i]+2);
			else if (argv[i][0] == 'Z')
				ePCR_default_pcr_size = atoi(argv[i]+2);
			else if (argv[i][0] == 'T')
			  ePCR_threads = atoi(argv[i]+2);
			else if (argv[i][0] == 'I')
				ePCR_iupac_mode = atoi(argv[i]+2);
			else if (argv[i][0] == 'X')
				three_prime_match = atoi(argv[i]+2);
		}
		else if (argv[i][0] == '-')    // -option
		{
			if (strcmp(argv[i],"-help") ==0)
				return Usage();
			if (strcmp(argv[i],"-margin") ==0)  // for backward compatibility
				margin = atoi(argv[++i]);  
		}
		else   // filename
		{
			if (stsfile == NULL)
				stsfile = argv[i];
			else if (seqfile ==NULL)
				seqfile = argv[i];
			else
				fprintf(stderr,"Argument \"%s\" ignored\n",argv[i]);
		}
	}

	if (ePCR_default_pcr_size < ePCR_DEFAULT_PCR_SIZE_MIN
	    || ePCR_default_pcr_size > ePCR_DEFAULT_PCR_SIZE_MAX) {
	  fprintf (stderr, "-Z argument is out of range (%u-%u)\n",
		   ePCR_DEFAULT_PCR_SIZE_MIN,  ePCR_DEFAULT_PCR_SIZE_MAX);
	  return Usage();
	}

	if (ePCR_quiet > 1 || ePCR_priority > 30 || ePCR_threads == 0 || ePCR_iupac_mode > 1) {
	  fprintf (stderr, "One or more of the Q=, I=, P=, or T= arguments is invalid\n");
	  return Usage();
	}
	
	if (strcasecmp(ePCR_outfile, "stdout") != 0) {
	  FILE *f = fopen (ePCR_outfile, "w");
	  if (!f) {
	    fprintf (stderr, "Error: can't open output file '%s'\n", ePCR_outfile);
	    stsfile = NULL;
	  }
	  fclose(f);	// Let ePCR_printf take care of the opens, when necessary
	}
	
	if (stsfile==NULL || seqfile==NULL)
		return Usage();
	
	///// Read STS primers database

	PCRmachine * e_PCR = new PCRmachine;

	e_PCR->SetWordSize(wdsize);
	e_PCR->SetMargin(margin);
	e_PCR->SetMismatch(mmatch);
	e_PCR->SetThreePrimeMatch(three_prime_match);

	if (!ePCR_quiet) {
		fprintf (stderr, "me-PCR parameters:\n");
		fprintf (stderr, "\twdsize=%d\n", e_PCR->GetWordSize());
		fprintf (stderr, "\tmargin=%d\n", e_PCR->GetMargin());
		fprintf (stderr, "\tmismatch=%d\n", e_PCR->GetMismatch());
		fprintf (stderr, "\tthree_prime_match=%d\n", e_PCR->GetThreePrimeMatch());
		fprintf (stderr, "\tquiet=%d (0=verbose mode; 1=no messages)\n", ePCR_quiet);
#ifdef __MWERKS__
		fprintf (stderr, "\tpriority=%d (range: %d (takes forever) to %d (locks machine))\n", 
						ePCR_priority, ePCR_PRIORITY_MIN, ePCR_PRIORITY_MAX);
#endif
		fprintf (stderr, "\toutfile=%s\n", ePCR_outfile);
		fprintf (stderr, "\tthreads=%d\n", ePCR_threads);
		fprintf (stderr, "\tmax STS line length=%d\n", ePCR_STS_line_length);
		fprintf (stderr, "\n");
	}

	
	if (ePCR_FileSize(seqfile) < e_PCR->MIN_FILESIZE_FOR_THREADING) {
	  if (!ePCR_quiet)
	    fprintf (stderr, "Notice: only one thread will be used because file is so small.\n");
	  ePCR_threads = 1;
	}

	if (!e_PCR->ReadStsFile(stsfile))
		return 1;

	///// Process sequence database (FASTA format)

	FastaFile fafile(SEQTYPE_NT);

	if (!fafile.Open(seqfile,"rb"))
		return 1;

	if (!ePCR_quiet)
	  fprintf (stderr, "m_margin=%d, max_pcr=%d\n", e_PCR->GetMargin(), e_PCR->max_pcr_size);

	// Fetch a vector of fasta seqs from fasta file
	fafile.Read();
	
	for (unsigned i=0; i<fafile.NumSeqs(); i++) {
	  FastaSeq **seqs = fafile.Seqs();
	  e_PCR->ProcessSeq(seqs[i]->Label(),seqs[i]->Sequence(), seqs[i]->Length());
	}

	fafile.Close();

#ifndef __MWERKS__
	/* The poor Mac SIOUX user _requires_ an indication that me-PCR has finished */
	if (!ePCR_quiet) {
#endif
	  if (ePCR_hits)
	    fprintf (stderr, "\n%lu %s found\n",
		     ePCR_hits, (ePCR_hits>1)?"hits":"hit");
	  else
	    fprintf (stderr, "\nNO HITS\n");
	  fprintf (stderr, "me-PCR complete; cleaning up ... \n\n###\n");
#ifndef __MWERKS__
	}
#endif

	delete e_PCR;  // not needed, but be tidy: close sts file, etc

	return 0;

} // end main
