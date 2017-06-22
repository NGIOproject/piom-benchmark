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

  int rank, size, dest, tag;
  int i, j;

  int nx, ny, nxp, nyp, xprocs, yprocs, barrier;

  MPI_Status status;

  int count, blocklength, stride;

  MPI_Datatype my_mpi_vector;

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

  /*
   *  Define the appropriate derived vector type
   *  The pattern is the same on all processes so everyone can make
   *  the same call
   */

  count = nxp;
  blocklength = nyp;
  stride = ny;

  MPI_Type_vector(count, blocklength, stride, MPI_FLOAT, &my_mpi_vector);

  /*  Commit it before use */

  MPI_Type_commit(&my_mpi_vector);


  if (rank == 0)
    {

      /*
       *  Send data to each of the workers in turn using synchronous sends
       *  Can do the sends in the usual order as we avoid the copy (except
       *  for the master)
       */

      for (dest = 1; dest < size; dest++)
	{

	  /*  Need to work out where to place the vector type within buf  */

	  istart = pcoords[dest][0]*nxp;
	  jstart = pcoords[dest][1]*nyp;

	  printf("rank %d sending to rank %d\n", rank, dest);

	  /*  Note that we send one vector */

	  MPI_Ssend(&(buf[istart][jstart]), 1, my_mpi_vector, dest, tag, comm);
	}

      /*  Master still needs to do a copy himself  */

      for (i=0; i < nxp; i++)
	{
	  for (j=0; j < nyp; j++)
	    {
	      x[i][j] = buf[i][j];
	    }
	}
    }
  else
    {
      /*
       *  Workers receive data from master
       */

      printf("rank %d receiving from rank 0\n", rank);

      /*  Note that we receive nxp*nyp reals */

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

  MPI_Type_free(&my_mpi_vector);

  MPI_Finalize();

  free(x);
  if(rank == 0){
    free(buf);
  }
  free(pcoords);
}
