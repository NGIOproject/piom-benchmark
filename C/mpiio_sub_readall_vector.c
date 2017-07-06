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
   *  x contains the local data only
   */
  float **x;

  int rank, size;
  int i, j;

  int nx, ny, nxp, nyp, xprocs, yprocs, subnum, ioworker, lowerrank, upperrank, barrier;

  int tag;

  MPI_File fh;
  MPI_Status status;
  MPI_Offset disp;

  int count, blocklength, stride, floatsize;

  MPI_Datatype my_mpi_vector;

  int istart, jstart;

  char filename[MAXFILENAME];

  double starttime, endtime, totaltime;

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm iocomm;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  checkandgetargumentssub(argc, argv, &nx, &ny, &xprocs, &yprocs, &nxp, &nyp, &subnum, &ioworker, &lowerrank, &upperrank, &barrier, size, rank);

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
   *  Initialise the arrays to a grey value
   */

  initarray(&x[0][0]  , nxp, nyp);

  if(barrier){
    MPI_Barrier(comm);
  }

  starttime = MPI_Wtime();

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


  /* Because file open is a collective operation, create a communicator for ioworkers and for non-ioworkers.*/
  /* This can then be used for a subset of all processes to do the I/O.*/
  MPI_Comm_split(comm, ioworker, rank, &iocomm);

  tag = 11;

  if(ioworker == 1){
    
    /*
     *  Open the file for reading only and attach to file handle fh
     *  No IO hints are passed since MPI_INFO_NULL is specified
     */
    
    if (MPI_File_open(iocomm, argv[1], MPI_MODE_RDONLY,
		      MPI_INFO_NULL, &fh) != MPI_SUCCESS)
      {
	printf("Open error on rank %d\n", rank);
      }

    MPI_Type_size(MPI_FLOAT, &floatsize);

    for(i=lowerrank+1;i<=upperrank;i++){

      istart = pcoords[i][0] * nxp;
      jstart = pcoords[i][1] * nyp;
      
      disp = istart;
      disp = disp*(MPI_Offset)ny;
      disp = disp + (MPI_Offset)jstart;
      disp = disp*(MPI_Offset)floatsize;

      /*
       *  Set view for this process using appropriate datatype
       */
      
      if (MPI_File_set_view(fh, disp, MPI_FLOAT, my_mpi_vector, "native",
			  MPI_INFO_NULL) != MPI_SUCCESS)
	{
	  printf("View error on rank %d\n", rank);
	}
      
      /*
       *  Read all the data for this process (ie nxp*nyp floats)
       */

      if (MPI_File_read_all(fh, &x[0][0], nxp*nyp, MPI_FLOAT, &status) != MPI_SUCCESS)
	{
	  printf("Read error on rank %d\n", rank);
	}


      MPI_Ssend(&x[0][0], nxp*nyp, MPI_FLOAT, i, tag, comm);
    }

    istart = pcoords[rank][0] * nxp;
    jstart = pcoords[rank][1] * nyp;

    disp = istart;
    disp = disp*(MPI_Offset)ny;
    disp = disp + (MPI_Offset)jstart;
    disp = disp*(MPI_Offset)floatsize;
            

    /*
     *  Set view for this process using appropriate datatype
     */
    
    if (MPI_File_set_view(fh, disp, MPI_FLOAT, my_mpi_vector, "native",
			  MPI_INFO_NULL) != MPI_SUCCESS)
      {
	printf("View error on rank %d\n", rank);
      }
    
      /*
       *  Read all the data for this process (ie nxp*nyp floats)
       */
    
    if (MPI_File_read_all(fh, &x[0][0], nxp*nyp, MPI_FLOAT, &status) != MPI_SUCCESS)
      {
	printf("Read error on rank %d\n", rank);
	}
    
    
    /*
     *  Close file
     */
    
    if (MPI_File_close(&fh) != MPI_SUCCESS){
	printf("Close error on rank %d\n", rank);
    }

  }else{
    // If I'm not an I/O worker receive our data
    MPI_Recv(&x[0][0], nxp*nyp, MPI_FLOAT, lowerrank, tag, comm, &status);    
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
  MPI_Comm_free(&iocomm);
  MPI_Finalize();
  free(x);
  free(pcoords);

}
