#piom-benchmark
The repository contains benchmark code for Parallel I/O models.  These are different ways of doing I/O from a parallel program (i.e. master I/O, individual files, collective MPI-I/O calls, etc...).

We have implemented the I/O models in both C and Fortran, to enable comparison of both programming languages.  We have also provided example submission scripts, designed to work on ARCHER (a Cray XC30 system).

This contains material contributed from a number of authors, including:

* Adrian Jackson, EPCC,
* David Henty, EPCC

## Structure
The repository has the following directory structure:
* C: The C source code
* F: The Fortran source code
* mkdata:  A code to generate input files for the benchmark



