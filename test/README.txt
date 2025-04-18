SIMPLE TESTING:
===============
Run 'me-PCR test.sts test.fa' to test me-PCR.  You should see one hit.

REGRESSION TESTING (for me-PCR developers):
===========================================

If you have Perl, you can run slow, exhaustive, and not very
enlightening tests by doing the following:

perl make_epcr_tests all
make test

Note that in the standard error output of the make, some trials will
generate the following error:

  WARNING: 1 STSs have a primer length sum greater than the pcr size: expected pcr size adjusted

This is benign; it's due to random adjustments made to the pcr size
given in the STS file.

Kevin Murphy
murphy@genome.chop.edu
Updated 2008 Jan 31

