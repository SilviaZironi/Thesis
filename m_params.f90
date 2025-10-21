module m_params
  use nrtype, only: dp
  implicit none
  private
  public :: get_params

contains

  subroutine get_params(v0, a, status, fname)
    ! legge i parametri che gli si passano da terminale, oppure usa quelli di params.in
    real(dp), intent(out) :: v0(2), a(2)
    integer,  intent(out), optional :: status
    character(*), intent(in),  optional :: fname

    integer :: narg, ios
    character(len=64)  :: s
    character(len=256) :: file

    narg = command_argument_count()

    if (narg >= 4) then
      call get_command_argument(1, s, status=ios)
      if (ios /= 0) call die('bad v01', status)
      read(s, *, iostat=ios) v0(1)
      if (ios /= 0) call die('bad v01', status)

      call get_command_argument(2, s, status=ios)
      if (ios /= 0) call die('bad v02', status)
      read(s, *, iostat=ios) v0(2)
      if (ios /= 0) call die('bad v02', status)

      call get_command_argument(3, s, status=ios)
      if (ios /= 0) call die('bad a1', status)
      read(s, *, iostat=ios) a(1)
      if (ios /= 0) call die('bad a1', status)

      call get_command_argument(4, s, status=ios)
      if (ios /= 0) call die('bad a2', status)
      read(s, *, iostat=ios) a(2)
      if (ios /= 0) call die('bad a2', status)

      if (present(status)) status = 0
      return
    end if

    if (present(fname)) then
      file = fname
    else
      file = 'params.in'
    end if

    call read_params_file(trim(file), v0, a, ios)
    if (ios /= 0) call die('failed reading '//trim(file)// &
                           ' (expected: v01 v02 a1 a2)', status)
    if (present(status)) status = 0
  contains
    subroutine die(msg, status)
      character(*), intent(in) :: msg
      integer, intent(out), optional :: status
      if (present(status)) then
        status = 1
      else
        write(*,*) 'FATAL: ', trim(msg)
        stop 2
      end if
    end subroutine
  end subroutine get_params


  subroutine read_params_file(fname, v0, a, ios)
    character(*), intent(in)  :: fname
    real(dp),     intent(out) :: v0(2), a(2)
    integer,      intent(out) :: ios
    integer :: u
    character(len=512) :: line

    ios = 0
    open(newunit=u, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) return

    do
      read(u,'(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      if (index(line,'!') > 0) line(index(line,'!'):)= ' '
      if (index(line,'#') > 0) line(index(line,'#'):)= ' '
      read(line,*, iostat=ios) v0(1), v0(2), a(1), a(2)
      if (ios == 0) exit
    end do
    close(u)
  end subroutine read_params_file

end module m_params
