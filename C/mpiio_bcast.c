#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

#include "ioutils.h"
#include "arralloc.h"

#define NDIM 2

int main(int argc, char **argv){
  /*
   *  pcoords stores the grid positions of each process
   */

  int **pcoords;
  
  /*
   *  buf is the large buffer for the master to read into
   *  x contains the local data only
   */

  float **buf;
  float **x;

  int rank, size;
  int i, j;

  int nx, ny, nxp, nyp, xprocs, yprocs, barrier;

  int istart, jstart;

  double starttime, endtime, totaltime;

  char filename[MAXFILENAME];

  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
  
  checkandgetarguments(argc, argv, &nx, &ny, &xprocs, &yprocs, &nxp, &nyp, &barrier, size, rank);

  pcoords = (int **)arralloc(sizeof(int), 2, size, NDIM);
  buf = (float **)arralloc(sizeof(float),2,nx,ny);
  x = (float **)arralloc(sizeof(float),2,nxp,nyp);

  /*
   *  Work out the coordinates of all the processes in the grid and
   *  print them out
   */
  initpgrid(&pcoords[0][0], xprocs, yprocs);

  if (rank == 0){
      printf("Running on %d process(es) in a %d x %d grid\n",
	     size, xprocs, yprocs);
      printf("\n");

      for (i=0; i < size; i++)
	{
	  printf("Process %2d has grid coordinates (%2d, %2d)\n",
		 i, pcoords[i][0], pcoords[i][1]);
	}
      printf("\n");
    }

  /*
   *  Initialise the arrays to a grey value
   */

  initarray(&buf[0][0], nx,  ny );
  initarray(&x[0][0]  , nxp, nyp);

  if(barrier){
    MPI_Barrier(comm);
  }

  starttime = MPI_Wtime();

  /*
   *  Read the entire array on the master process
   */

  if (rank == 0){
      ioread (argv[1], &buf[0][0], nx*ny);
      printf("\n");
    }

  /*
   *  Broadcast the data
   */

  MPI_Bcast(&buf[0][0], nx*ny, MPI_FLOAT, 0, comm);

  endtime = MPI_Wtime();

  /*
   *  Copy down the correct data from buf to x. Need to work out, using the
   *  position in the process grid, what the index of the bottom-left-hand
   *  pixel for this process.
   */

  istart = pcoords[rank][0]*nxp;
  jstart = pcoords[rank][1]*nyp;

  for (i=0; i < nxp; i++){
      for (j=0; j < nyp; j++){
	  x[i][j] = buf[istart+i][jstart+j];
	}
    }

  /*
   *  Every process writes out its local data array x to an individually
   *  named file which has the rank appended to the file name
   */

#ifdef DEBUG
  createfilename(filename, "coutput", nxp, nyp, rank);
  iowrite(filename, &x[0][0], nxp*nyp);
#endif

  totaltime = endtime - starttime;

  dotimings(totaltime, rank, size);

  MPI_Finalize();

  free(x);
  free(buf);
  free(pcoords);
}
