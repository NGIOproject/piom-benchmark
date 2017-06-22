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

  int istart, jstart;

  double starttime, endtime, totaltime;

  MPI_Status status;

  int sizes[NDIM], subsizes[NDIM], starts[NDIM];

  MPI_Datatype *my_mpi_subarray;

  char filename[MAXFILENAME];

  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  checkandgetarguments(argc, argv, &nx, &ny, &xprocs, &yprocs, &nxp, &nyp, &barrier, size, rank);

  my_mpi_subarray = malloc(size*sizeof(MPI_Datatype));
  pcoords = (int **)arralloc(sizeof(int), 2, size, NDIM);
  if(rank == 0){
    buf = (float **)arralloc(sizeof(float),2,nx,ny);
  }
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
   *  Define the correct nxp x nyp subarray for each of the processes
   *  in the calculation
   */

  sizes[0] = nx;
  sizes[1] = ny;

  subsizes[0] = nxp;
  subsizes[1] = nyp;

  for (dest = 0; dest < size; dest++){
    starts[0] = pcoords[dest][0] * nxp;
    starts[1] = pcoords[dest][1] * nyp;
    
    MPI_Type_create_subarray(NDIM, sizes, subsizes, starts,
			     MPI_ORDER_C, MPI_FLOAT, &my_mpi_subarray[dest]);
    
      MPI_Type_commit(&my_mpi_subarray[dest]);
  }
  
  if (rank == 0){

    /*
     *  Send data to each of the workers in turn using synchronous sends
     *  Can do the sends in the usual order as we avoid the copy (except
     *  for the master)
     */
    
    for (dest = 1; dest < size; dest++)
      {
	printf("rank %d sending to rank %d\n", rank, dest);
	
	/*
	 *  Note that we send one subarray
	 *  We pass the starting address of buf as the position of the subarray
	 *  within buf is already encoded in the deifinition of my_mpi_subarray
	 */
	
	MPI_Ssend(&x[0][0], 1, my_mpi_subarray[dest], dest, tag, comm);
      }
    
    /*  Master still needs to do a copy himself  */
    
    for (i=0; i < nxp; i++)
      {
	for (j=0; j < nyp; j++)
	  {
	    x[i][j] = buf[i][j];
	  }
      }
  }else{
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

  for (dest = 0; dest < size; dest++){

      MPI_Type_free(&my_mpi_subarray[dest]);
  }

  MPI_Finalize();
  free(my_mpi_subarray);
  free(x);
  if(rank == 0){
    free(buf);
  }
  free(pcoords);

}
