#include <stdio.h>
#include <stdlib.h>
#include "ioutils.h"
#include "mpi.h"


void checkandgetargumentssub(int argc, char **argv, int *nx, int *ny, int *xprocs, int *yprocs, int *nxp, int *nyp, int *subnum, int *ioworker, int *lowerrank, int *upperrank, int *barrier, int size, int rank){
  
  int chunk;

  if ( argc != 8 ){
    /* We print argv[0] assuming it is the program name */
    if(rank == 0){
      printf("usage: %s filename nx ny xprocs yprocs barrier numio", argv[0] );
      printf("This application expects you to provide the input file now, the size of the input data set (nx*ny), the number of processes you want in each dimension (xproc and yproc), and whether to use a barrier before timing (barrier, 0 is no barrier, 1 uses a barrier), the number of I/O workers you want (numio)\n");
    }
    MPI_Finalize();
    exit(-1);
    
  }else{
    getargs(argc,argv,nx,ny,xprocs,yprocs,nxp,nyp,barrier,size,rank);
    *subnum = atoi(argv[7]);
    if((size/(*subnum))*(*subnum) != size){
      if(rank == 0){
	printf("The number of I/O workers specified (%d) does not exactly divide the number of processes being used (%d)", *subnum, size);
      }
      MPI_Finalize();
      exit(-1);      
    }
    // Work out how many processes each I/O worker is reponsible for
    chunk = size/(*subnum);
    // Work out if this process is an I/O worker
    if(rank%chunk == 0){
      *ioworker = 1;
      // If we are an I/O worker then work out the range of ranks we are responsible for (including itself)
      *lowerrank = rank;
      *upperrank = rank + chunk - 1;
    
    }else{
      *ioworker = 0;
      // If we are not an I/O worker then record which process is our responsible I/O worker
      *lowerrank = rank-(rank%chunk);
      *upperrank = rank-(rank%chunk);
    }

  }
  if(rank == 0){
     printf("number of i/o workers: %d\n",*subnum);
  }
}

void checkandgetarguments(int argc, char **argv, int *nx, int *ny, int *xprocs, int *yprocs, int *nxp, int *nyp, int *barrier, int size, int rank){

  if ( argc != 7 ){
    /* We print argv[0] assuming it is the program name */
    if(rank == 0){
      printf("usage: %s filename nx ny xprocs yprocs barrier", argv[0] );
      printf("This application expects you to provide the input file now, the size of the input data set (nx*ny), the number of processes you want in each dimension (xproc and yproc), and whether to use a barrier before timing (barrier, 0 is no barrier, 1 uses a barrier)\n");
    }
    MPI_Finalize();
    exit(-1);
    
  }else{
    getargs(argc,argv,nx,ny,xprocs,yprocs,nxp,nyp,barrier,size,rank);
  }

}
void getargs(int argc, char **argv, int *nx, int *ny, int *xprocs, int *yprocs, int *nxp, int *nyp, int *barrier, int size, int rank){

  *nx = atoi(argv[2]);
  *ny = atoi(argv[3]);
  *xprocs = atoi(argv[4]);
  *yprocs = atoi(argv[5]);
  *barrier = atoi(argv[6]);
  
  if((*xprocs)*(*yprocs) != size){
    if(rank == 0){
      printf("Number of processes being used (%d) does not match number requested in the two dimension (%d x %d)\n",size,*xprocs,*yprocs);
    }
    MPI_Finalize();
    exit(-1);
  }
  
  *nxp = (*nx)/(*xprocs);
  *nyp = (*ny)/(*yprocs);
  
  if((*nxp)*(*xprocs) != *nx){
    if(rank == 0){
      printf("Number of processes in x dimension (xprocs) does not exactly divide the x dimension of the input (nx).  Quitting\n");
    }
    MPI_Finalize();
    exit(-1);
  }
  
  if((*nyp)*(*yprocs) != *ny){
    if(rank == 0){
      printf("Number of processes in y dimension (yprocs) does not exactly divide the y dimension of the input (ny).  Quitting\n");
    }
    MPI_Finalize();
    exit(-1);
  }
  
  if(*barrier != 0 && *barrier != 1){
    if(rank == 0){
      printf("Barrier should be 0 or 1\n");
    }
    MPI_Finalize();
    exit(-1);
  }

  if(rank == 0){
    printf("Running on %d processes\n",size);
    printf("nx: %d ny: %d nxp: %d nyp: %d xprocs: %d yprocs: %d\n",*nx,*nx,*nxp,*nyp,*xprocs,*yprocs);
    printf("barrier: %d\n",*barrier);
  }

}

void dotimings(double totaltime, int rank, int size){
  
  double avtotaltime, mintotaltime, maxtotaltime;

  MPI_Reduce(&totaltime, &avtotaltime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  avtotaltime = avtotaltime/size;

  MPI_Reduce(&totaltime, &maxtotaltime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&totaltime, &mintotaltime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

  if(rank == 0){
    printf("Times (seconds): Average: %lf Maximum: %lf Minimum: %lf\n",avtotaltime,maxtotaltime,mintotaltime);
  }


}

void createfilename(char *filename, char *basename, int nx, int ny, int rank)
{
  if (rank < 0)
    {
      sprintf(filename, "%s%04dx%04d.dat", basename, nx, ny);
    }
  else
    {
      sprintf(filename, "%s%04dx%04d_%02d.dat", basename, nx, ny, rank);
    }
}


void iosize(char *filename, int *nx, int *ny)
{ 
  char *tmpname;

  tmpname = filename;

  while (*tmpname < '0' || *tmpname > '9')
  {
    tmpname++;
  }

  if (sscanf(tmpname, "%04dx%04d.dat", nx, ny) != 2)
  {
    printf("iosize: error parsing filename <%s>\n", filename);
    exit(1);
  }
}


void ioread(char *filename, void *ptr, int nfloat)
{
  int i;

  FILE *fp;
  
  float *data = (float *) ptr;

  printf("ioread: reading <%s> ...\n", filename);

  if ( (fp = fopen(filename, "r")) == NULL)
  {
    printf("ioread: failed to open input file <%s>\n", filename);
    exit(1);
  }

  if (fread(data, sizeof(float), nfloat, fp) != nfloat) 
  {
    printf("ioread: error reading input file <%s>\n", filename);
    exit(1);
  }

  fclose(fp);

  printf("... done\n", filename);
}

void iochunkread(char *filename, void *ptr, int nfloat, long int offset)
{
  
  FILE *fp;
  
  float *data = (float *) ptr;

//  printf("iochunkread: reading <%s> ...\n", filename);

  if ( (fp = fopen(filename, "r")) == NULL)
  {
    printf("iochunkread: failed to open input file <%s>\n", filename);
    exit(1);
  }

  /* Initial datsa offset (moving to the correct starting point in the file */
  if (fseek(fp, offset, SEEK_SET) != 0){
    printf("iochunkread: failed to move the the correct place in the input file <%s>\n", filename);
    exit(1);
  }
  
  if (fread(data, sizeof(float), nfloat, fp) != nfloat) 
    {
      printf("iochunkread: error reading input file <%s>\n", filename);
      exit(1);
    }

  fclose(fp);

//  printf("... done\n", filename);
}

void iowrite(char *filename, void *ptr, int nfloat)
{
  int i;

  FILE *fp;
  
  float *data = (float *) ptr;

  printf("iowrite: writing <%s> ...\n", filename);

  if ( (fp = fopen(filename, "w")) == NULL)
  {
    printf("iowrite: failed to open output file <%s>\n", filename);
    exit(1);
  }

  if (fwrite(data, sizeof(float), nfloat, fp) != nfloat) 
  {
    printf("iowrite: error writing output file <%s>\n", filename);
    exit(1);
  }

  fclose(fp);
  printf("... done\n", filename);
}



void iochunkwrite(char *filename, void *ptr, int nfloat, long int offset)
{
  int i;

  FILE *fp;
  
  float *data = (float *) ptr;

  printf("iochunkwrite: writing <%s> ...\n", filename);

  if ( (fp = fopen(filename, "w")) == NULL)
  {
    printf("iochunkwrite: failed to open output file <%s>\n", filename);
    exit(1);
  }

  if (fseek(fp, offset, SEEK_SET) != 0){
    printf("iochunkread: failed to move the the correct place in the input file <%s>\n", filename);
    exit(1);
  }

  if (fwrite(data, sizeof(float), nfloat, fp) != nfloat) 
  {
    printf("iowchunkrite: error writing output file <%s>\n", filename);
    exit(1);
  }

  fclose(fp);
  printf("... done\n", filename);
}


#define INITDATAVAL 0.5

void initarray(void *ptr, int nx, int ny)
{
  int i, j;

  float *data = (float *) ptr;

  for (i=0; i < nx*ny; i++)
    {
      data[i] = INITDATAVAL;
    }
}

#define NDIM 2

void initpgrid(void *ptr, int nxproc, int nyproc)
{
  MPI_Comm gridcomm;
  MPI_Comm comm = MPI_COMM_WORLD;

  int dims[NDIM];
  int periods[NDIM] = {0, 0};

  int reorder = 0;

  int *pcoords = (int *) ptr;

  int i;

  dims[0] = nxproc;
  dims[1] = nyproc;

  MPI_Cart_create(comm, NDIM, dims, periods, reorder, &gridcomm);

  for (i=0; i < nxproc*nyproc; i++)
    {
      MPI_Cart_coords(gridcomm, i, NDIM, pcoords+NDIM*i);
    }

  MPI_Comm_free(&gridcomm);

}
