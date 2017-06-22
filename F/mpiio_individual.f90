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
!  x contains the local data only
!

  real, dimension(:,:), allocatable :: x

  integer :: rank, size, ierr
  integer :: i, j

  integer :: istart, jstart

  character*(maxfilename) :: filename

  double precision :: starttime, endtime, totaltime

  integer :: comm = MPI_COMM_WORLD
  integer(kind=8) :: offset

  call MPI_INIT(ierr)

  call MPI_COMM_SIZE(comm, size, ierr)
  call MPI_COMM_RANK(comm, rank, ierr)

  call checkandgetarguments(filename, maxfilename, nx, ny, xprocs, yprocs, nxp, nyp, size, rank)

  allocate(pcoords(ndim,size))
  allocate(x(nxp,nyp))
  
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

  call initarray(x,   nxp, nyp)

  starttime = mpi_wtime()

!
!  Read the entire array on the master process
!  Passing "-1" as the rank argument means that the file name has no
!  trailing "_rank" appended to it, ie we read the global file
!


  offset = pcoords(2,rank+1)*nx
  offset = offset*nyp
  offset = offset + pcoords(1,rank+1)*nyp
  do i=1,nyp
    call iochunkread (filename, x(:,i), nxp, offset)
    offset = offset + nx
 end do

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

  deallocate(pcoords,x)

end

