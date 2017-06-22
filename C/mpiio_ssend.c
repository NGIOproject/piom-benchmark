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
   *  buf is the large buffer for the master to read into
   *  x contains the local data only
   */
  float **buf;
  float **x;

  int rank, size, dest, tag;
  int i, j;

  int nx, ny, nxp, nyp, xprocs, yprocs, barrier;

  MPI_Status status;

  int istart, jstart;

  double starttime, endtime, totaltime;

  char filename[MAXFILENAME];

  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  checkandgetarguments(argc, argv, &nx, &ny, &xprocs, &yprocs, &nxp, &nyp, &barrier, size, rank);

  pcoords = (int **)arralloc(sizeof(int), 2, size, NDIM);
  if(rank == 0){
    buf = (float **)arralloc(sizeof(float),2,nx,ny);
  }
  x = (float **)arralloc(sizeof(float),2,nxp,nyp);

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
   *  Initialise the arrays to a grey value
   */

  if(rank == 0){
    initarray(&buf[0][0], nx,  ny );
  }
  initarray(&x[0][0]  , nxp, nyp);

  if(barrier){
    MPI_Barrier(comm);
  }

  starttime = MPI_Wtime();

  /*
   *  Read the entire array on the master process
   *  Passing "-1" as the rank argument means that the file name has no
   *  trailing "_rank" appended to it, ie we read the global file
   */

  if (rank == 0)
    {
      ioread (argv[1], &buf[0][0], nx*ny);
      printf("\n");
    }

  tag = 0;

  if (rank == 0)
    {

      /*
       *  Send data to each of the workers in turn using synchronous sends
       *  Do them in reverse order of rank so that the code is slightly
       *  neater - the master does the copy but omits the send to itself.
       *  Tag all messages with tag=0
       */

      for (dest = size-1; dest >= 0; dest--)
	{
	  

  /*
   *  Copy down the correct data from buf to x. Need to work out, using the
   *  position in the process grid, what the index of the bottom-left-hand
   *  pixel for this particular destiniation process.
   */

	  istart = pcoords[dest][0]*nxp;
	  jstart = pcoords[dest][1]*nyp;

	  for (i=0; i < nxp; i++)
	    {
	      for (j=0; j < nyp; j++)
		{
		  x[i][j] = buf[istart+i][jstart+j];
		}
	    }

	  if (dest != 0)
	    {
	      printf("rank %d sending to rank %d\n", rank, dest);
	      MPI_Ssend(&x[0][0], nxp*nyp, MPI_FLOAT, dest, tag, comm);
	    }
	}
    }
  else
    {
      /*
       *  Workers receive data from master
       */

      printf("rank %d receiving from rank 0\n", rank);
      MPI_Recv(&x[0][0], nxp*nyp, MPI_FLOAT, 0, tag, comm, &status);
    }

  endtime = MPI_Wtime();

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
  if(rank == 0){
    free(buf);
  }
  free(pcoords);

}
