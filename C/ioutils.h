/*
 *  The maximum length of a file name
 */

#define MAXFILENAME 200

void initarray(void *ptr, int nx, int ny);
void initpgrid(void *ptr, int nxproc, int nyproc);

void checkandgetargumentssub(int argc, char **argv, int *nx, int *ny, int *xprocs, int *yprocs, int *nxp, int *nyp, int *subnum, int *ioworker, int *lowerrank, int *upperrank, int *barrier, int size, int rank);
void checkandgetarguments(int argc, char **argv, int *nx, int *ny, int *xprocs, int *yprocs, int *nxp, int *nyp, int *barrier, int size, int rank);
void getargs(int argc, char **argv, int *nx, int *ny, int *xprocs, int *yprocs, int *nxp, int *nyp, int *barrier, int size, int rank);

void dotimings(double totaltime, int rank,int size);

void createfilename(char *filename, char *basename, int nx, int ny, int rank);

void ioread (char *filename, void *ptr, int nfloat);
void iowrite(char *filename, void *ptr, int nfloat);

void iochunkread (char *filename, void *ptr, int nfloat, long int offset);
void iochunkwrite(char *filename, void *ptr, int nfloat, long int offset);
