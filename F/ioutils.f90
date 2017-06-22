subroutine createfilename(filename, basename, nx, ny, rank)

  implicit none

  character*(*) :: filename, basename

  integer :: nx, ny, rank

  if (rank .lt. 0) then

    write(filename, fmt='(i10.10,''x'',i10.10,''.dat'')') nx, ny

  else

    write(filename, fmt='(i10.10,''x'',i10.10,''_'',i4.4,''.dat'')')  nx, ny, rank

  end if

  filename = basename//filename

end

subroutine iosize(filename, nx, ny)

  implicit none

  character*(*) :: filename
  integer       :: nx, ny

  integer :: i

  do i = 1, len(filename)

    if (iachar(filename(i:i)) .ge. iachar('0') .and. &
        iachar(filename(i:i)) .le. iachar('9')         ) exit

  end do

  if (i > len(filename)) then

    write(*,*) 'iosize: error parsing filename ', filename

    nx = -1
    ny = -1

  else

    read(filename(i:len(filename)), fmt='(i4.4,''x'',i4.4,''.dat'')') nx, ny

  end if

end



subroutine ioread(filename, data, nreal)

  implicit none

  integer, parameter :: iounit   = 10
  integer, parameter :: realsize = 4

  character*(*) :: filename
  integer       :: nreal
  real          :: data(nreal)

  integer :: i

  write(*,*) 'ioread: reading ', filename

  open(unit=iounit, file=filename, form='unformatted', &
       access='direct', recl=realsize)

  do i = 1, nreal
    read(unit=iounit, rec=i) data(i)
  end do

  close(unit=iounit)

  write(*,*) '.. done'

end

subroutine iochunkread(filename, data, nreal, offset)

  implicit none

  integer, parameter :: iounit   = 10
  integer, parameter :: realsize = 4

  character*(*) :: filename
  integer       :: nreal
  real          :: data(nreal)
  integer(kind=8)::offset

  integer(kind=8) :: i

!AJ remove large amounts of output for large files
!AJ  write(*,*) 'iochunkread: reading ', filename

  open(unit=iounit, file=filename, form='unformatted', &
       access='direct', recl=realsize)

  do i = 1,nreal     
    read(unit=iounit, rec=i+offset) data(i)
  end do

  close(unit=iounit)

!AJ remove large amount of output for large files
!AJ  write(*,*) '.. done'

end


subroutine iowrite(filename, data, nreal)

  implicit none

  integer, parameter :: iounit   = 10
  integer, parameter :: realsize = 4

  character*(*) :: filename
  integer       :: nreal
  real          :: data(nreal)

  integer :: i

  write(*,*) 'iowrite: writing ', filename

  open(unit=iounit, file=filename, form='unformatted', &
       access='direct', recl=realsize)

  do i = 1, nreal
    write(unit=iounit, rec=i) data(i)
  end do

  close(unit=iounit)

  write(*,*) '.. done'

end

subroutine initarray(data, nx, ny)

  implicit none

  real, parameter :: initdataval = 0.5

  integer :: nx, ny

  real data(nx*ny)

  integer :: i

  do i = 1, nx*ny
    data(i) = initdataval
  end do

end

subroutine initpgrid(pcoords, nxproc, nyproc)

  use mpi

  implicit none

  integer, parameter :: ndim = 2

  integer :: nxproc, nyproc
  integer, dimension(ndim, nxproc*nyproc) :: pcoords

  integer, dimension(ndim) :: dims, periods

  integer :: i, ierr
  integer :: comm = MPI_COMM_WORLD
  integer :: gridcomm

  logical :: reorder

  periods = (/ 0, 0 /)
  reorder = .false.

  dims(1) = nxproc
  dims(2) = nyproc

  call MPI_CART_CREATE(comm, ndim, dims, periods, reorder, gridcomm, ierr)

  do i = 1, nxproc*nyproc

    call MPI_CART_COORDS(gridcomm, i-1, ndim, pcoords(1, i), ierr)

  end do

  call MPI_COMM_FREE(gridcomm, ierr)

end

subroutine checkandgetarguments(filename, maxfilename, nx, ny, xprocs, yprocs, nxp, nyp, size, rank)

  use mpi

  implicit none

  integer, intent(in) :: size, rank, maxfilename
  integer, intent(out) :: nx, ny, xprocs, yprocs, nxp, nyp
  character*(maxfilename), intent(out) :: filename

  integer :: i, numargs, ierr
  character(len=32) :: arg
 
  numargs = command_argument_count()
  if(numargs .ne. 5) then
     if(rank .eq. 0) then
        write(*,*) "usage: inputfilename nx ny xprocs yproc"
        write(*,*) "This application expects you to provide the input file now, the size of the input data set (nx*ny) and the number of processes you want in each dimension (xproc and yproc)"
     end if
     call MPI_FINALIZE(ierr)
     stop
  end if
  call getarg(1,filename)
  call getarg(2,arg)
  arg = trim(arg)
  read(arg, '(I10)') nx
  call getarg(3,arg)
  arg = trim(arg)
  read(arg, '(I10)') ny
  call getarg(4,arg)
  arg = trim(arg)
  read(arg, '(I10)') xprocs
  call getarg(5,arg)
  arg = trim(arg)
  read(arg, '(I10)') yprocs
  if(xprocs*yprocs .ne. size) then
     if(rank .eq. 0) then
        write(*,*) 'The specified xprocs and yprocs assignment does not match the total number of processes being used.'
        write(*,*) 'xprocs is ',xprocs,' yprocs is ',yprocs,' but total processes used is ',size,' which does not match xprocs * yprocs'
     end if
     call MPI_FINALIZE(ierr)
     stop
  end if
  nxp = nx/xprocs
  if(xprocs*nxp .ne. nx) then
     if(rank .eq. 0) then
        write(*,*) 'nx does not exactly divide by xprocs.  Stopping the code!'
     end if
     call MPI_FINALIZE(ierr)
     stop
  end if
  nyp = ny/yprocs
  if(yprocs*nyp .ne. ny) then
     if(rank .eq. 0) then
        write(*,*) 'ny does not exactly divide by yprocs. Stopping the code!'
     end if
     call MPI_FINALIZE(ierr)
     stop
  end if
  
end subroutine checkandgetarguments


subroutine dotimings(totaltime, rank, size)

  use mpi

  implicit none

  double precision, intent(in) :: totaltime
  integer, intent(in) :: rank, size
  double precision :: avtotaltime, mintotaltime, maxtotaltime
  integer :: ierr

  call MPI_Reduce(totaltime, avtotaltime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  avtotaltime = avtotaltime/size

  call MPI_Reduce(totaltime, maxtotaltime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(totaltime, mintotaltime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD, ierr)

  if(rank .eq. 0) then
   write(*,*) "Times (seconds): Average: ",avtotaltime," Maximum: ",maxtotaltime," Minimum: ",mintotaltime
  end if

end subroutine dotimings
