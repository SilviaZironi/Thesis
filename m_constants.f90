module constants
use nrtype
implicit none
public

! Input/output constants

  integer(kind=i1b), parameter :: inp=1
  integer(kind=i1b), parameter :: oup=2

! Mathematical constants

  real(kind=dp), parameter :: half=0.50_dp
  real(kind=dp), parameter :: zero=0.0_dp
  real(kind=dp), parameter :: one=1.0_dp
  real(kind=dp), parameter :: two=2.0_dp
  real(kind=dp), parameter :: three=3.0_dp
  real(kind=dp), parameter :: four=4.0_dp
  real(kind=dp), parameter :: five=5.0_dp
  real(kind=dp), parameter :: six=6.0_dp
  real(kind=dp), parameter :: seven=7.0_dp
  real(kind=dp), parameter :: eight=8.0_dp
  real(kind=dp), parameter :: nine=9.0_dp
  real(kind=dp), parameter :: ten=10.0_dp
  complex(kind=dp), parameter :: imat=(zero,one)

! Physical constants

  real(kind=dp), parameter :: hbarc=197.3269788_dp
  real(kind=dp), parameter :: m=938.91897_dp
  real(kind=dp), parameter :: mf=m/hbarc
  real(kind=dp), parameter :: de=-0.1_dp

end module constants