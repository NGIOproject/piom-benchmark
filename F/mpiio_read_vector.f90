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

  real, dimension(:,:), allocatable :: x

  integer :: rank, size, ierr
  integer :: i, j

  character*(maxfilename) :: filename

  double precision :: starttime, endtime, totaltime

!
!  Variables needed for MPI-IO
!

  integer, dimension(MPI_STATUS_SIZE) :: status
  integer :: count, blocklength, stride, istart, jstart, barrier
  integer :: realsize
  integer :: my_mpi_vector, fh
  integer (kind=MPI_OFFSET_KIND) :: disp = 0

  integer :: comm = MPI_COMM_WORLD

  call MPI_INIT(ierr)

  call MPI_COMM_SIZE(comm, size, ierr)
  call MPI_COMM_RANK(comm, rank, ierr)

  call checkandgetarguments(filename, maxfilename, nx, ny, xprocs, yprocs, nxp, nyp, barrier, size, rank)

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

  if(barrier .eq. 1) then
     call MPI_Barrier(comm, ierr)
  end if

  starttime = mpi_wtime()

!
!  Define the NXP x NYP vector for this distribution
!  Note that it is the same for every process. To ensure that each
!  process reads different data from the file, they use different disps
!

  count = nyp
  blocklength = nxp
  stride = nx

  call MPI_TYPE_VECTOR(count, blocklength, stride,    &
                       MPI_REAL, my_mpi_vector, ierr)

!  Commit it before use

  call MPI_TYPE_COMMIT(my_mpi_vector, ierr)

!
!  Open the file for reading only and attach to file handle fh
!  No IO hints are passed since MPI_INFO_NULL is specified
!

  call MPI_FILE_OPEN(comm, filename, MPI_MODE_RDONLY, &
                     MPI_INFO_NULL, fh, ierr)

  if (ierr /= MPI_SUCCESS) write(*,*) 'Open error on rank ', rank

!
!  Set view for this process using the vector type with an appropriate
!  value if disp (computed in BYTES!)
!

  istart = pcoords(1, rank+1)*nxp
  jstart = pcoords(2, rank+1)*nyp

  call MPI_TYPE_SIZE(MPI_REAL, realsize, ierr)

  disp = (jstart*nx + istart)*realsize

  call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL, my_mpi_vector, 'native', &
                         MPI_INFO_NULL, ierr)

  if (ierr /= MPI_SUCCESS) write(*,*) 'View error on rank ', rank

!
!  Read all the data for this process (ie NXP*NYP floats)
!

  call MPI_FILE_READ(fh, x, nxp*nyp, MPI_REAL, status, ierr)

  if (ierr /= MPI_SUCCESS) write(*,*) 'Read error on rank ', rank

!
!  Close file
!

  call MPI_FILE_CLOSE(fh, ierr)

  if (ierr /= MPI_SUCCESS) write(*,*) 'Close error on rank ', rank

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

