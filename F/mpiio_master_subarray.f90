program mpiio

  use mpi

  implicit none

!
!  The global data size is nx x ny
!

  integer :: nx
  integer :: ny

!
!  The processes are in a 2D array of dimension XPROCS x YPROCS, with
!  a total of NPROCS processes
!

  integer, parameter :: ndim = 2

  integer :: xprocs
  integer :: yprocs

!
!  The local data size is NXP x NYP
!

  integer :: nxp
  integer :: nyp

!
!  The maximum length of a file name
!

  integer, parameter :: maxfilename = 64

!
!  pcoords stores the grid positions of each process
!

  integer, dimension(:,:), allocatable :: pcoords

!
!  buf is the large buffer for the master to read into
!  x contains the local data only
!

  real, dimension(:,:), allocatable :: buf
  real, dimension(:,:), allocatable :: x

  integer :: rank, size, ierr, dest, tag, barrier
  integer :: i, j

  integer, dimension(MPI_STATUS_SIZE) :: status

!  Variables needed for subarray types

  integer, dimension(ndim)   :: sizes, subsizes, starts
  integer, dimension(:), allocatable :: my_mpi_subarray

  character*(maxfilename) :: filename

  double precision :: starttime, endtime, totaltime

  integer :: comm = MPI_COMM_WORLD

  call MPI_INIT(ierr)

  call MPI_COMM_SIZE(comm, size, ierr)
  call MPI_COMM_RANK(comm, rank, ierr)

  call checkandgetarguments(filename, maxfilename, nx, ny, xprocs, yprocs, nxp, nyp, barrier, size, rank)

  allocate(pcoords(ndim,size))
  if(rank .eq. 0) then
    allocate(buf(nx,ny))
  end if
  allocate(x(nxp,nyp))
  allocate(my_mpi_subarray(size))
!
!  Work out the coordinates of all the processes in the grid and
!  print them out
!

  call initpgrid(pcoords, xprocs, yprocs)

  if (rank .eq. 0) then

     write(*,*) 'Running on ', size, ' process(es) in a ', &
                xprocs, ' x ', yprocs, ' grid'
     write(*,*)

     do i = 0, size-1
       write(*,*) 'Process ', i, ' has grid coordinates (', &
                  pcoords(1,i+1), ', ', pcoords(2,i+1), ')'
     end do

     write(*,*)

  end if
  
!
!  Initialise the arrays to a grey value
!
  if(rank .eq. 0) then
    call initarray(buf, nx,  ny )
  end if
  call initarray(x,   nxp, nyp)

  if(barrier .eq. 1) then
     call MPI_Barrier(comm, ierr)
  end if

  starttime = mpi_wtime()
!
!  Read the entire array on the master process
!  Passing "-1" as the rank argument means that the file name has no
!  trailing "_rank" appended to it, ie we read the global file
!

  if (rank .eq. 0) then

    call ioread (filename, buf, nx*ny)
    write(*,*)

  end if

  tag = 0

!
!  Define the correct NXP x NYP subarray for each of the processes
!  in the calculation
!

  sizes(1) = nx
  sizes(2) = ny

  subsizes(1) = nxp
  subsizes(2) = nyp

  do dest = 0, size-1

     starts(1) = pcoords(1, dest+1) * nxp
     starts(2) = pcoords(2, dest+1) * nyp;

     call MPI_TYPE_CREATE_SUBARRAY(ndim, sizes, subsizes, starts,  &
                                   MPI_ORDER_FORTRAN, MPI_REAL, &
                                   my_mpi_subarray(dest+1), ierr)
!  Commit it before use


     call MPI_TYPE_COMMIT(my_mpi_subarray(dest+1), ierr)

  end do

  if (rank == 0) then 

!  Send data to each of the workers in turn using synchronous sends
!  Can do the sends in the usual order as we avoid the copy (except
!  for the master)

     do dest = 1, size-1

        write(*,*) 'rank ', rank, ' sending to rank ', dest

!  Note that we send one subarray
!  We pass the starting address of buf as the position of the subarray
!  within buf is already encoded in the deifinition of my_mpi_subarray

        call MPI_SSEND(buf, 1, my_mpi_subarray(dest+1), &
                       dest, tag, comm, ierr)

     end do

!  Master still needs to do a copy to himself
     do j = 1, nyp
        do i = 1, nxp

           x(i,j) = buf(i,j)
              
        end do
     end do

  else

!  Workers receive data from master

     write(*,*) 'rank ', rank, ' receiving from rank 0'

!  Note that we receive nxp*nyp reals

     call MPI_RECV(x, nxp*nyp, MPI_REAL, 0, tag, comm, status, ierr)

  end if

  endtime = mpi_wtime()

!
!  Every process writes out its local data array x to an individually
!  named file which as the rank appended to the file name
!

!  call createfilename(filename, 'foutput', nxp, nyp, rank)
!  call iowrite(filename, x, nxp*nyp)

  totaltime = endtime-starttime
  call dotimings(totaltime, rank, size)

  call MPI_FINALIZE(ierr)

  deallocate(pcoords,x,buf,my_mpi_subarray)
  if(rank .eq. 0) then
    deallocate(buf)
  end if

end

