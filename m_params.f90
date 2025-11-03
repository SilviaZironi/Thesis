module m_params
  use nrtype, only: dp
  implicit none
  private
  public :: get_params

contains

  subroutine get_params(v0, a, status)
    real(dp), intent(out) :: v0(2), a(2)
    integer,  intent(out), optional :: status
    integer :: narg, ios

    v0 = 0.0_dp; a = 0.0_dp
    narg = command_argument_count()
    if (narg < 4) then
      call die('Usage: trimod v01 v02 a1 a2', status); return
    end if

    call parse_arg(1,'v01',v0(1),ios); if (ios/=0) then; call die('bad v01',status); return; end if
    call parse_arg(2,'v02',v0(2),ios); if (ios/=0) then; call die('bad v02',status); return; end if
    call parse_arg(3,'a1', a(1),ios);  if (ios/=0) then; call die('bad a1', status); return; end if
    call parse_arg(4,'a2', a(2),ios);  if (ios/=0) then; call die('bad a2', status); return; end if

    if (present(status)) status = 0

  contains
    subroutine parse_arg(i, name, val, ios)
      integer,  intent(in)  :: i
      character(*), intent(in) :: name
      real(dp), intent(out) :: val
      integer,  intent(out) :: ios
      character(len=64) :: s
      call get_command_argument(i, s, status=ios); if (ios/=0) return
      call normalize_decimal(s)
      read(s, *, iostat=ios) val
    end subroutine

    subroutine normalize_decimal(s)
      character(*), intent(inout) :: s
      integer :: k
      do k = 1, len_trim(s)
        if (s(k:k)==',') s(k:k)='.'
      end do
    end subroutine

    subroutine die(msg, status)
      character(*), intent(in) :: msg
      integer, intent(out), optional :: status
      if (present(status)) then
        status = 1
      else
        write(*,*) trim(msg); stop 1
      end if
    end subroutine
  end subroutine get_params

end module m_params


