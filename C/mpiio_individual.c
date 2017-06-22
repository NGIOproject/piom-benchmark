#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

#include "ioutils.h"

#include "arralloc.h"

#define NDIM 2

int main(int argc, char **argv)
{
  /*
   *  pcoords stores the grid positions of each process
   */

  int **pcoords;
  
  /*
   *  x contains the local data only
   */
  float **x;

  int rank, size;
  int i, j;

  int nx, ny, nxp, nyp, xprocs, yprocs, barrier;

  MPI_Status status;

  int count, blocklength, stride;

  MPI_Datatype my_mpi_vector;

  int istart, jstart;

  double starttime, endtime, totaltime;

  char filename[MAXFILENAME];

  MPI_Comm comm = MPI_COMM_WORLD;

  long int offset;
  int datasize;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  checkandgetarguments(argc, argv, &nx, &ny, &xprocs, &yprocs, &nxp, &nyp, &barrier, size, rank);


  pcoords = (int **)arralloc(sizeof(int), 2, size, NDIM);
  x = (float **)arralloc(sizeof(float),2,nxp,nyp);

  /*
   *  Work out the coordinates of all the processes in the grid and
   *  print them out
   */

  initpgrid(&pcoords[0][0], xprocs, yprocs);

  if (rank == 0)
    {
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
   *  Initialise the array to a grey value
   */
  initarray(&x[0][0]  , nxp, nyp);


  /*
   *  Read the entire array on the master process
   *  Passing "-1" as the rank argument means that the file name has no
   *  trailing "_rank" appended to it, ie we read the global file
   */

  if(barrier){
    MPI_Barrier(comm);
  }

  starttime = MPI_Wtime();

  datasize = sizeof(float);

  offset = pcoords[rank][0];
  offset = offset*ny;
  offset = offset*nxp;
  offset = offset + pcoords[rank][1]*nyp;
  offset = offset*datasize;
  for(i=0; i<nxp; i++){
    iochunkread (argv[1], &x[i][0], nyp, offset);
    offset = offset + ny*datasize;
  }

  printf("\n");

  endtime = MPI_Wtime();

#ifdef DEBUG
  createfilename(filename, "coutput", nxp, nyp, rank);
  iowrite(filename, &x[0][0], nxp*nyp);
#endif

  totaltime = endtime - starttime;

  dotimings(totaltime, rank, size);

  MPI_Finalize();

  free(x);
  free(pcoords);
}
