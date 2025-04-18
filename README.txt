me-PCR version 1.0.X

INSTALLING BINARIES
===================

*If* you see a 'bin' subfolder, then this package of me-PCR contains
executable binaries for a few platforms.  Installation is currently
manual, by copying me-PCR (or, for Windows, me-PCR.exe and
pthreadVC.dll) into a directory in your search path.  Alternatively,
you can run me-PCR from any directory by changing into that directory.
For UNIXish operations systems, you may need to invoke it as
"./me-PCR" instead of just "me-PCR".  For Windows users, more
information is available in the README file in the windows subfolder
in the bin folder.

COMPILING
=========
To compile and install me-PCR, you need two things:
1) a C++ compiler
2) a POSIX threads library

For UNIXish operating systems, just type 'cd src; make' to compile the
program.  You will see some warnings because the code has not been
cleaned up yet.  

There are special makefiles for AIX and the Borland 5.5 compiler under
Windows: makefile.aix and makefile.borland, respectively.  To use
these, e.g., run 'cd src; make -f makefile.borland'.

If you encounter any errors, please contact mepcr@genome.chop.edu.
Our lab runs me-PCR on AIX and Mac OS X, so there may be problems with
some platforms.  We're happy to work through these problems with you
and improve the package.  Eventually we will have a GNU configure
rig set up.

Installation is manual; just copy the program where you would like it.
A manual is provided in -mdoc (man) and HTML formats.

Various platform-specific notes follow:

-------
WINDOWS
-------
To compile on Windows, you have two options at the moment, both free:

1) Use the freely available Borland C++ compiler.  Go to
http://www.borland.com/products/downloads/download_cbuilder.html
and click on the 'Compiler' link.  You will need to establish an account.
Use 'make -f makefile.borland' to compile me-PCR and link 
with the threads DLL that we include in the borland subdirectory.  You will need to
copy the DLL pthreadVC.dll to the same directory as the generated exectutable.

2) Download and install Cygwin from http://www.cygwin.com/.  You will
need to make sure to download the C++ (GNU C++) compiler and GNU make.
me-PCR should compile fine with the standard makefile.

--------
AIX
--------

For AIX, you need to add '-Xlinker -bmaxdata:0x80000000' to the e-PCR link command
(the command for the e-PCR target).

For 64-bit AIX kernel with 32-bit GNU tools (hmm), you need to export
OBJECT_MODE=32 before compiling.


###
2003-05-02 Kevin Murphy <murphy@genome.chop.edu>
Updated 2008-02-18 Kevin Murphy
