module trimod_sub
implicit none
private
public :: gauss,potential1,potential2,channels,permutation,tsingle,tcoupled,search

  interface gauss
          module procedure gauss
  end interface

  interface potential1
          module procedure potential1
  end interface

  interface potential2
          module procedure potential2
  end interface

  interface channels
          module procedure channels
  end interface

  interface permutation
          module procedure permutation
  end interface

  interface tsingle
          module procedure tsingle
  end interface

  interface tcoupled
          module procedure tcoupled
  end interface

  interface search
          module procedure search
  end interface

contains
!================================================================================================================================!
subroutine gauss(a,b,c,n,x,x2,wx)
use nrtype
use nr_sub
implicit none
!--------------------------------------------------------------------!
! This subroutine determines Gaussian quadrature points for the
! interval [a,c] where half of the points lie to the left and
! the right of b, respectively.
! The subroutine gauleg from Numerical Recipes determines the n
! gauss points gp and the weights gw in the interval [-1.0,1.0]
!--------------------------------------------------------------------!
integer(kind=i4b), intent(in) :: n
real(kind=dp), intent(in) :: a,b,c
real(kind=dp), dimension(n), intent(out) :: x,x2,wx
integer(kind=i4b) :: i
real(kind=dp) :: b1,b2,b3,b4,b5,b6,cn
real(kind=dp), dimension(n) :: gp,gw

  call gauleg(-1.0_dp,1.0_dp,gp,gw)
  b1=2.0_dp*a*c-a*b-b*c
  b2=2.0_dp*a*(b-c)
  b3=a-2.0_dp*b+c
  b4=a-c
  b5=b-c
  b6=b-a
  do i=1,n
     cn=b3*gp(i)+b4
     x(i)=(b1*(1.0_dp+gp(i))+b2)/cn
     x2(i)=x(i)*x(i)
     wx(i)=2.0_dp*b5*b4*b6*gw(i)/cn**2
  end do

end subroutine gauss
!================================================================================================================================!
subroutine potential1(tta,jja,interaction,ngh,k,pot)
use nrtype
use constants
use cdbonn_potential
use chn4lo500_2017
implicit none
!---------------------------------------------------------------------------------------------!
! Set up the potential in momentum space, units fm^2
!---------------------------------------------------------------------------------------------!
integer(kind=i4b), intent(in) :: tta,jja,interaction,ngh
real(kind=dp), dimension(ngh), intent(in) :: k
real(kind=dp), dimension(6,ngh,ngh), intent(out) :: pot
integer(kind=i4b) :: i,j
real(kind=dp) :: qk,qkp
real(kind=dp), dimension(6) :: vv

! Generate the NN potential

select case (interaction)
   case (1)
      do i=1,ngh
         qk = k(i)*hbarc
         if (qk == zero) then
            qk=0.000001_dp
         end if
         do j=1,ngh
            qkp=k(j)*hbarc
            if (qkp == zero) then
               qkp=0.000001_dp
            end if
            call potcdbonn(jja,tta,qk,qkp,vv)
            pot(1,i,j) = vv(1)*hbarc**2
            pot(2,i,j) = vv(2)*hbarc**2
            pot(3,i,j) = vv(3)*hbarc**2
            pot(4,i,j) = vv(4)*hbarc**2
            pot(5,i,j) = vv(5)*hbarc**2
            pot(6,i,j) = vv(6)*hbarc**2
         end do
      end do
   case (2)
      do i=1,ngh
         qk = k(i)*hbarc
         if (qk == zero) then
            qk=0.000001_dp
         end if
         do j=1,ngh
            qkp=k(j)*hbarc
            if (qkp == zero) then
               qkp=0.000001_dp
            end if
            call potn4lo500_2017(jja,tta,qk,qkp,vv)
            pot(1,i,j) = vv(1)*hbarc**2
            pot(2,i,j) = vv(2)*hbarc**2
            pot(3,i,j) = vv(3)*hbarc**2
            pot(4,i,j) = vv(4)*hbarc**2
            pot(5,i,j) = vv(5)*hbarc**2
            pot(6,i,j) = vv(6)*hbarc**2
         end do
      end do
   case default
      write(*,*) "Bad input for interaction in potential1"
      stop
end select

end subroutine potential1
!================================================================================================================================!
subroutine potential2(force,jmom,ostat,cutnumr,ngh,k,pot)
use nrtype
use constants
use ichiral
implicit none
!---------------------------------------------------------------------------------------------!
!  Set up the potential in k space, units MeV**(-2)
!
!  The INPUT parameters are OSTAT, FORCE, PPUN1, PPUN2, JMOM and CUTNUMR
!
!  FORCE= 'nn', 'np', 'pp'
!  OSTAT=0:  LO
!  OSTAT=1:  NLO
!  OSTAT=2:  NNLO
!  OSTAT=3:  NNNLO
!  OSTAT=4:  NNNNLO
!  PPUN1:    final center-of-mass momentum of the nucleons in GeV
!  PPUN2:    initial center-of-mass momentum of the nucleons in GeV
!  CUTNUMR:  cut-off used in the LS equation
!  JMOM:     total angular momentum
!
!  The output vector is vv and contain the NN potential.
!  In order to bring the PWD MEs of this routine into the convention of "chiral_n3lo.f90",
!  multiply vv(i) by the factor of PI/2D0/1000000D0. 
!
!  vv(1):    S=0, L1=L2=J
!  vv(2):    S=1, L1=L2=J
!  vv(3):    S=1, L1=L2=J-1
!  vv(4):    S=1, L1=J-1, L2=J+1
!  vv(5):    S=1, L1=J+1, L2=J-1
!  vv(6):    S=1, L1=L2=J+1
!---------------------------------------------------------------------------------------------!
  character(len=2), intent(in) :: force
  integer(kind=i4b), intent(in) :: jmom,ostat,cutnumr,ngh
  real(kind=dp), dimension(ngh), intent(in) :: k
  real(kind=dp), dimension(6,ngh,ngh), intent(out) :: pot
  real(kind=dp) :: ppun1,ppun2
  real(kind=dp), dimension(6) :: vv
  integer(kind=i4b) :: i,j

  do i=1,ngh
     ppun1 = k(i)*hbarc/(ten**3)
     if (ppun1 == zero) then
        ppun1=0.000001_dp
     end if
     do j=1,ngh
        ppun2=k(j)*hbarc/(ten**3)
        if (ppun2 == zero) then
           ppun2=0.000001_dp
        end if
        call chiralmompwd(ostat,force,ppun1,ppun2,jmom,cutnumr,vv)
        pot(1,i,j) = vv(1)*(hbarc**2)/(ten**6)
        pot(2,i,j) = vv(2)*(hbarc**2)/(ten**6)
        pot(3,i,j) = vv(6)*(hbarc**2)/(ten**6)
        pot(4,i,j) = vv(3)*(hbarc**2)/(ten**6)
        pot(5,i,j) = vv(5)*(hbarc**2)/(ten**6)
        pot(6,i,j) = vv(4)*(hbarc**2)/(ten**6)
     end do
  end do

end subroutine potential2
!================================================================================================================================!
subroutine channels(nchmax,chan)
use nrtype
implicit none
integer(kind=i4b), intent(in) :: nchmax
integer(kind=i4b), dimension(nchmax,6), intent(out) :: chan
!--------------------------------------------------------------------!
! This subroutine sets up the channels
! The quantum numbers are all doubled
! The convention is: ( l , s , j , t , lambda , I )
!--------------------------------------------------------------------!

! Channel 1   ( 0 , 0 , 0 , 1 , 0 , 1/2 )   1S0
  chan(1,1)=0 ; chan(1,2)=0 ; chan(1,3)=0 ; chan(1,4)=2 ; chan(1,5)=0 ; chan(1,6)=1

! Channel 2   ( 1 , 1 , 0 , 1 , 1 , 1/2 )   3P0
  chan(2,1)=2 ; chan(2,2)=2 ; chan(2,3)=0 ; chan(2,4)=2 ; chan(2,5)=2 ; chan(2,6)=1

! Channel 3   ( 1 , 0 , 1 , 0 , 1 , 1/2 )   1P1
  chan(3,1)=2 ; chan(3,2)=0 ; chan(3,3)=2 ; chan(3,4)=0 ; chan(3,5)=2 ; chan(3,6)=1

! Channel 4   ( 1 , 0 , 1 , 0 , 1 , 3/2 )   1P1
  chan(4,1)=2 ; chan(4,2)=0 ; chan(4,3)=2 ; chan(4,4)=0 ; chan(4,5)=2 ; chan(4,6)=3

! Channel 5   ( 1 , 1 , 1 , 1 , 1 , 1/2 )   3P1
  chan(5,1)=2 ; chan(5,2)=2 ; chan(5,3)=2 ; chan(5,4)=2 ; chan(5,5)=2 ; chan(5,6)=1

! Channel 6   ( 1 , 1 , 1 , 1 , 1 , 3/2 )   3P1
  chan(6,1)=2 ; chan(6,2)=2 ; chan(6,3)=2 ; chan(6,4)=2 ; chan(6,5)=2 ; chan(6,6)=3

! Channel 7   ( 0 , 1 , 1 , 0 , 0 , 1/2 )   3S1
  chan(7,1)=0 ; chan(7,2)=2 ; chan(7,3)=2 ; chan(7,4)=0 ; chan(7,5)=0 ; chan(7,6)=1

! Channel 8   ( 2 , 1 , 1 , 0 , 0 , 1/2 )   3D1
  chan(8,1)=4 ; chan(8,2)=2 ; chan(8,3)=2 ; chan(8,4)=0 ; chan(8,5)=0 ; chan(8,6)=1

! Channel 9   ( 0 , 1 , 1 , 0 , 2 , 3/2 )   3S1
  chan(9,1)=0 ; chan(9,2)=2 ; chan(9,3)=2 ; chan(9,4)=0 ; chan(9,5)=4 ; chan(9,6)=3

! Channel 10  ( 2 , 1 , 1 , 0 , 2 , 3/2 )   3D1
  chan(10,1)=4 ; chan(10,2)=2 ; chan(10,3)=2 ; chan(10,4)=0 ; chan(10,5)=4 ; chan(10,6)=3

! Channel 11  ( 2 , 0 , 2 , 1 , 2 , 3/2 )   1D2
  chan(11,1)=4 ; chan(11,2)=0 ; chan(11,3)=4 ; chan(11,4)=2 ; chan(11,5)=4 ; chan(11,6)=3

! Channel 12  ( 2 , 0 , 2 , 1 , 2 , 5/2 )   1D2
  chan(12,1)=4 ; chan(12,2)=0 ; chan(12,3)=4 ; chan(12,4)=2 ; chan(12,5)=4 ; chan(12,6)=5

! Channel 13  ( 2 , 1 , 2 , 0 , 2 , 3/2 )   3D2
  chan(13,1)=4 ; chan(13,2)=2 ; chan(13,3)=4 ; chan(13,4)=0 ; chan(13,5)=4 ; chan(13,6)=3

! Channel 14  ( 2 , 1 , 2 , 0 , 2 , 5/2 )   3D2
  chan(14,1)=4 ; chan(14,2)=2 ; chan(14,3)=4 ; chan(14,4)=0 ; chan(14,5)=4 ; chan(14,6)=5

! Channel 15  ( 1 , 1 , 2 , 1 , 1 , 3/2 )   3P2
  chan(15,1)=2 ; chan(15,2)=2 ; chan(15,3)=4 ; chan(15,4)=2 ; chan(15,5)=2 ; chan(15,6)=3

! Channel 16  ( 3 , 1 , 2 , 1 , 1 , 3/2 )   3F2
  chan(16,1)=6 ; chan(16,2)=2 ; chan(16,3)=4 ; chan(16,4)=2 ; chan(16,5)=2 ; chan(16,6)=3

! Channel 17  ( 1 , 1 , 2 , 1 , 3 , 5/2 )   3P2
  chan(17,1)=2 ; chan(17,2)=2 ; chan(17,3)=4 ; chan(17,4)=2 ; chan(17,5)=6 ; chan(17,6)=5

! Channel 18  ( 3 , 1 , 2 , 1 , 3 , 5/2 )   3F2
  chan(18,1)=6 ; chan(18,2)=2 ; chan(18,3)=4 ; chan(18,4)=2 ; chan(18,5)=6 ; chan(18,6)=5

! Channel 19  ( 3 , 0 , 3 , 0 , 3 , 5/2 )   1F3
  chan(19,1)=6 ; chan(19,2)=0 ; chan(19,3)=6 ; chan(19,4)=0 ; chan(19,5)=6 ; chan(19,6)=5

! Channel 20  ( 3 , 0 , 3 , 0 , 3 , 7/2 )   1F3
  chan(20,1)=6 ; chan(20,2)=0 ; chan(20,3)=6 ; chan(20,4)=0 ; chan(20,5)=6 ; chan(20,6)=7

! Channel 21  ( 3 , 1 , 3 , 1 , 3 , 5/2 )   3F3
  chan(21,1)=6 ; chan(21,2)=2 ; chan(21,3)=6 ; chan(21,4)=2 ; chan(21,5)=6 ; chan(21,6)=5

! Channel 22  ( 3 , 1 , 3 , 1 , 3 , 7/2 )   3F3
  chan(22,1)=6 ; chan(22,2)=2 ; chan(22,3)=6 ; chan(22,4)=2 ; chan(22,5)=6 ; chan(22,6)=7

! Channel 23  ( 2 , 1 , 3 , 0 , 2 , 5/2 )   3D3
  chan(23,1)=4 ; chan(23,2)=2 ; chan(23,3)=6 ; chan(23,4)=0 ; chan(23,5)=4 ; chan(23,6)=5

! Channel 24  ( 4 , 1 , 3 , 0 , 2 , 5/2 )   3G3
  chan(24,1)=8 ; chan(24,2)=2 ; chan(24,3)=6 ; chan(24,4)=0 ; chan(24,5)=4 ; chan(24,6)=5

! Channel 25  ( 2 , 1 , 3 , 0 , 4 , 7/2 )   3D3
  chan(25,1)=4 ; chan(25,2)=2 ; chan(25,3)=6 ; chan(25,4)=0 ; chan(25,5)=8 ; chan(25,6)=7

! Channel 26  ( 4 , 1 , 3 , 0 , 4 , 7/2 )   3G3
  chan(26,1)=8 ; chan(26,2)=2 ; chan(26,3)=6 ; chan(26,4)=0 ; chan(26,5)=8 ; chan(26,6)=7

! Channel 27  ( 4 , 0 , 4 , 1 , 4 , 7/2 )   1G4
  chan(27,1)=8 ; chan(27,2)=0 ; chan(27,3)=8 ; chan(27,4)=2 ; chan(27,5)=8 ; chan(27,6)=7

! Channel 28  ( 4 , 0 , 4 , 1 , 4 , 9/2 )   1G4
  chan(28,1)=8 ; chan(28,2)=0 ; chan(28,3)=8 ; chan(28,4)=2 ; chan(28,5)=8 ; chan(28,6)=9

! Channel 29  ( 4 , 1 , 4 , 0 , 4 , 7/2 )   3G4
  chan(29,1)=8 ; chan(29,2)=2 ; chan(29,3)=8 ; chan(29,4)=0 ; chan(29,5)=8 ; chan(29,6)=7

! Channel 30  ( 4 , 1 , 4 , 0 , 4 , 9/2 )   3G4
  chan(30,1)=8 ; chan(30,2)=2 ; chan(30,3)=8 ; chan(30,4)=0 ; chan(30,5)=8 ; chan(30,6)=9

! Channel 31  ( 3 , 1 , 4 , 1 , 3 , 7/2 )   3F4
  chan(31,1)=6 ; chan(31,2)=2 ; chan(31,3)=8 ; chan(31,4)=2 ; chan(31,5)=6 ; chan(31,6)=7

! Channel 32  ( 5 , 1 , 4 , 1 , 3 , 7/2 )   3H4
  chan(32,1)=10 ; chan(32,2)=2 ; chan(32,3)=8 ; chan(32,4)=2 ; chan(32,5)=6 ; chan(32,6)=7

! Channel 33  ( 3 , 1 , 4 , 1 , 5 , 9/2 )   3F4
  chan(33,1)=6 ; chan(33,2)=2 ; chan(33,3)=8 ; chan(33,4)=2 ; chan(33,5)=10 ; chan(33,6)=9

! Channel 34  ( 5 , 1 , 4 , 1 , 5 , 9/2 )   3H4
  chan(34,1)=10 ; chan(34,2)=2 ; chan(34,3)=8 ; chan(34,4)=2 ; chan(34,5)=10 ; chan(34,6)=9

! Channel 35  ( 5 , 0 , 5 , 0 , 5 , 9/2 )   1H5
  chan(35,1)=10 ; chan(35,2)=0 ; chan(35,3)=10 ; chan(35,4)=0 ; chan(35,5)=10 ; chan(35,6)=9

! Channel 36  ( 5 , 0 , 5 , 0 , 5 , 11/2 )   1H5
  chan(36,1)=10 ; chan(36,2)=0 ; chan(36,3)=10 ; chan(36,4)=0 ; chan(36,5)=10 ; chan(36,6)=11

! Channel 37  ( 5 , 1 , 5 , 1 , 5 , 9/2 )   3H5
  chan(37,1)=10 ; chan(37,2)=2 ; chan(37,3)=10 ; chan(37,4)=2 ; chan(37,5)=10 ; chan(37,6)=9

! Channel 38  ( 5 , 1 , 5 , 1 , 5 , 11/2 )   3H5
  chan(38,1)=10 ; chan(38,2)=2 ; chan(38,3)=10 ; chan(38,4)=2 ; chan(38,5)=10 ; chan(38,6)=11

! Channel 39  ( 4 , 1 , 5 , 0 , 4 , 9/2 )   3G5
  chan(39,1)=8 ; chan(39,2)=2 ; chan(39,3)=10 ; chan(39,4)=0 ; chan(39,5)=8 ; chan(39,6)=9

! Channel 40  ( 6 , 1 , 5 , 0 , 4 , 9/2 )   3I5
  chan(40,1)=12 ; chan(40,2)=2 ; chan(40,3)=10 ; chan(40,4)=0 ; chan(40,5)=8 ; chan(40,6)=9

! Channel 41  ( 4 , 1 , 5 , 0 , 6 , 11/2 )   3G5
  chan(41,1)=8 ; chan(41,2)=2 ; chan(41,3)=10 ; chan(41,4)=0 ; chan(41,5)=12 ; chan(41,6)=11

! Channel 42  ( 6 , 1 , 5 , 0 , 6 , 11/2 )   3I5
  chan(42,1)=12 ; chan(42,2)=2 ; chan(42,3)=10 ; chan(42,4)=0 ; chan(42,5)=12 ; chan(42,6)=11

! Channel 43  ( 6 , 0 , 6 , 1 , 6 , 11/2 )   1I6
  chan(43,1)=12 ; chan(43,2)=0 ; chan(43,3)=12 ; chan(43,4)=2 ; chan(43,5)=12 ; chan(43,6)=11

! Channel 44  ( 6 , 0 , 6 , 1 , 6 , 13/2 )   1I6
  chan(44,1)=12 ; chan(44,2)=0 ; chan(44,3)=12 ; chan(44,4)=2 ; chan(44,5)=12 ; chan(44,6)=13

! Channel 45  ( 6 , 1 , 6 , 0 , 6 , 11/2 )   3I6
  chan(45,1)=12 ; chan(45,2)=2 ; chan(45,3)=12 ; chan(45,4)=0 ; chan(45,5)=12 ; chan(45,6)=11

! Channel 46  ( 6 , 1 , 6 , 0 , 6 , 13/2 )   3I6
  chan(46,1)=12 ; chan(46,2)=2 ; chan(46,3)=12 ; chan(46,4)=0 ; chan(46,5)=12 ; chan(46,6)=13

! Channel 47  ( 5 , 1 , 6 , 1 , 5 , 11/2 )   3H6
  chan(47,1)=10 ; chan(47,2)=2 ; chan(47,3)=12 ; chan(47,4)=2 ; chan(47,5)=10 ; chan(47,6)=11

! Channel 48  ( 7 , 1 , 6 , 1 , 5 , 11/2 )   3J6
  chan(48,1)=14 ; chan(48,2)=2 ; chan(48,3)=12 ; chan(48,4)=2 ; chan(48,5)=10 ; chan(48,6)=11

! Channel 49  ( 5 , 1 , 6 , 1 , 7 , 13/2 )   3H6
  chan(49,1)=10 ; chan(49,2)=2 ; chan(49,3)=12 ; chan(49,4)=2 ; chan(49,5)=14 ; chan(49,6)=13

! Channel 50  ( 7 , 1 , 6 , 1 , 7 , 13/2 )   3J6
  chan(50,1)=14 ; chan(50,2)=2 ; chan(50,3)=12 ; chan(50,4)=2 ; chan(50,5)=14 ; chan(50,6)=13

! Channel 51  ( 7 , 0 , 7 , 0 , 7 , 13/2 )   1J7
  chan(51,1)=14 ; chan(51,2)=0 ; chan(51,3)=14 ; chan(51,4)=0 ; chan(51,5)=14 ; chan(51,6)=13

! Channel 52  ( 7 , 0 , 7 , 0 , 7 , 15/2 )   1J7
  chan(52,1)=14 ; chan(52,2)=0 ; chan(52,3)=14 ; chan(52,4)=0 ; chan(52,5)=14 ; chan(52,6)=15

! Channel 53  ( 7 , 1 , 7 , 1 , 7 , 13/2 )   3J7
  chan(53,1)=14 ; chan(53,2)=2 ; chan(53,3)=14 ; chan(53,4)=2 ; chan(53,5)=14 ; chan(53,6)=13

! Channel 54  ( 7 , 1 , 7 , 1 , 7 , 15/2 )   3J7
  chan(54,1)=14 ; chan(54,2)=2 ; chan(54,3)=14 ; chan(54,4)=2 ; chan(54,5)=14 ; chan(54,6)=15

! Channel 55  ( 6 , 1 , 7 , 0 , 6 , 13/2 )   3I7
  chan(55,1)=12 ; chan(55,2)=2 ; chan(55,3)=14 ; chan(55,4)=0 ; chan(55,5)=12 ; chan(55,6)=13

! Channel 56  ( 8 , 1 , 7 , 0 , 6 , 13/2 )   3K7
  chan(56,1)=16 ; chan(56,2)=2 ; chan(56,3)=14 ; chan(56,4)=0 ; chan(56,5)=12 ; chan(56,6)=13

! Channel 57  ( 6 , 1 , 7 , 0 , 8 , 15/2 )   3I7
  chan(57,1)=12 ; chan(57,2)=2 ; chan(57,3)=14 ; chan(57,4)=0 ; chan(57,5)=16 ; chan(57,6)=15

! Channel 58  ( 8 , 1 , 7 , 0 , 8 , 15/2 )   3K7
  chan(58,1)=16 ; chan(58,2)=2 ; chan(58,3)=14 ; chan(58,4)=0 ; chan(58,5)=16 ; chan(58,6)=15

! Channel 59  ( 8 , 0 , 8 , 1 , 8 , 15/2 )   1K8
  chan(59,1)=16 ; chan(59,2)=0 ; chan(59,3)=16 ; chan(59,4)=2 ; chan(59,5)=16 ; chan(59,6)=15

! Channel 60  ( 8 , 0 , 8 , 1 , 8 , 17/2 )   1K8
  chan(60,1)=16 ; chan(60,2)=0 ; chan(60,3)=16 ; chan(60,4)=2 ; chan(60,5)=16 ; chan(60,6)=17

! Channel 61  ( 8 , 1 , 8 , 0 , 8 , 15/2 )   3K8
  chan(61,1)=16 ; chan(61,2)=2 ; chan(61,3)=16 ; chan(61,4)=0 ; chan(61,5)=16 ; chan(61,6)=15

! Channel 62  ( 8 , 1 , 8 , 0 , 8 , 17/2 )   3K8
  chan(62,1)=16 ; chan(62,2)=2 ; chan(62,3)=16 ; chan(62,4)=0 ; chan(62,5)=16 ; chan(62,6)=17

! Channel 63  ( 7 , 1 , 8 , 1 , 7 , 15/2 )   3J8
  chan(63,1)=14 ; chan(63,2)=2 ; chan(63,3)=16 ; chan(63,4)=2 ; chan(63,5)=14 ; chan(63,6)=15

! Channel 64  ( 9 , 1 , 8 , 1 , 7 , 15/2 )   3L8
  chan(64,1)=18 ; chan(64,2)=2 ; chan(64,3)=16 ; chan(64,4)=2 ; chan(64,5)=14 ; chan(64,6)=15

! Channel 65  ( 7 , 1 , 8 , 1 , 9 , 17/2 )   3J8
  chan(65,1)=14 ; chan(65,2)=2 ; chan(65,3)=16 ; chan(65,4)=2 ; chan(65,5)=18 ; chan(65,6)=17

! Channel 66  ( 9 , 1 , 8 , 1 , 9 , 17/2 )   3L8
  chan(66,1)=18 ; chan(66,2)=2 ; chan(66,3)=16 ; chan(66,4)=2 ; chan(66,5)=18 ; chan(66,6)=17

end subroutine channels
!================================================================================================================================!
subroutine permutation(nchmax,nq,nx,chan,q,x,gfac)
use nrtype
use constants
use nr_sub
implicit none
integer(kind=i4b), parameter :: ikmax=24
integer(kind=i4b), intent(in) :: nchmax,nq,nx
integer(kind=i4b), dimension(nchmax,6), intent(in) :: chan
real(kind=dp), dimension(nq), intent(in) :: q
real(kind=dp), dimension(nx), intent(in) :: x
real(kind=dp), dimension(nchmax,nchmax,nq,nq,nx), intent(out) :: gfac
integer(kind=i4b) :: a,ap,iq,iqp,ix,ik,il1,ilp1,l,lp,l2,lp2
integer(kind=i4b) :: llmin,llmax,fmin,fmax,fpmin,fpmax,ll,ss,f,fp
real(kind=dp) :: sumk,prodf,prods,prodl
real(kind=dp) :: gg,xg,sixj,aninej
real(kind=dp), dimension(0:ikmax,nx) :: leg
real(kind=dp), dimension(nx) :: sumff,sumls
real(kind=dp), dimension(nq,nq,nx) :: sumq

  call factorial

! Precalculate Legendre polynomials

  do ix=1,nx
     do ik=0,ikmax
        leg(ik,ix) = plgndr_s(ik,0,x(ix))
     end do
  end do

! Calculate G_{a,ap} (q1,q2,x)

  do a=1,nchmax
     l=chan(a,1)/2
     do ap=1,nchmax
        lp=chan(ap,1)/2

        sumq=zero
        do il1=0,l
           do ilp1=0,lp

              l2=l-il1
              lp2=lp-ilp1

              gg = - sqrt(real(chan(a,1),kind=dp)+one)  * sqrt(real(chan(a,2),kind=dp)+one) &
                   * sqrt(real(chan(a,3),kind=dp)+one)  * sqrt(real(chan(a,4),kind=dp)+one) &
                   * sqrt(real(chan(a,5),kind=dp)+one)  * sqrt(real(chan(a,6),kind=dp)+one) &
                   * sqrt(real(chan(ap,1),kind=dp)+one) * sqrt(real(chan(ap,2),kind=dp)+one) &
                   * sqrt(real(chan(ap,3),kind=dp)+one) * sqrt(real(chan(ap,4),kind=dp)+one) &
                   * sqrt(real(chan(ap,5),kind=dp)+one) * sqrt(real(chan(ap,6),kind=dp)+one) &
                   * sqrt(gamma(real(chan(a,1),kind=dp)+two)/(gamma(two*il1+one)*gamma(two*l2+one))) &
                   * sqrt(gamma(real(chan(ap,1),kind=dp)+two)/(gamma(two*ilp1+one)*gamma(two*lp2+one))) &
                   * sixj(1,1,chan(a,4),1,1,chan(ap,4)) * half**(l2+ilp1)

              llmin=max(abs(chan(a,1)/2-chan(a,5)/2),abs(chan(ap,1)/2-chan(ap,5)/2))
              llmax=min(chan(a,1)/2+chan(a,5)/2,chan(ap,1)/2+chan(ap,5)/2)
              sumls=zero
              do ss=1,3,2
                 prods = (real(ss,kind=dp)+one) * sixj(1,1,chan(a,2),1,ss,chan(ap,2))
                          
                 do ll=llmin,llmax
                    prodl = (two*real(ll,kind=dp)+one) &
                          * aninej(chan(a,1),chan(a,2),chan(a,3),chan(a,5),1,chan(a,6),2*ll,ss,1) &
                          * aninej(chan(ap,1),chan(ap,2),chan(ap,3),chan(ap,5),1,chan(ap,6),2*ll,ss,1)

                    fmin=min(abs(l2-chan(a,5)/2),abs(il1-ll))
                    fmax=max(l2+chan(a,5)/2,il1+ll)
                    sumff=zero
                    do f=fmin,fmax
                       fpmin = min(abs(ilp1-chan(ap,5)/2),abs(lp2-ll))
                       fpmax = max(ilp1+chan(ap,5)/2,lp2+ll)
                       prodf = xg(2*l2,0,chan(a,5),0,2*f,0) * sixj(2*il1,2*l2,chan(a,1),chan(a,5),2*ll,2*f)
                       do fp=fpmin,fpmax

                          do ix=1,nx
                             sumk=zero
                             do ik=0,2*max(l,lp)+2
                                sumk = sumk + (two*real(ik,kind=dp)+one) &
                                     * xg(2*ik,0,2*lp2,0,2*f,0) * xg(2*ik,0,2*il1,0,2*fp,0) &
                                     * sixj(2*f,2*il1,2*ll,2*fp,2*lp2,2*ik) * leg(ik,ix) !* plgndr_s(ik,0,x(ix))
                             end do
                             sumff(ix) = sumff(ix) + prodf * sumk * xg(2*ilp1,0,chan(ap,5),0,2*fp,0) &
                                   * sixj(2*lp2,2*ilp1,chan(ap,1),chan(ap,5),2*ll,2*fp)
                          end do

                       end do
                    end do

                    do ix=1,nx
                       sumls(ix) = sumls(ix) + prods * prodl * sumff(ix)
                    end do

                 end do
              end do

              do iq=1,nq
                 do iqp=1,nq
                    do ix=1,nx
                       sumq(iq,iqp,ix) = sumq(iq,iqp,ix) + (q(iq)**(l-il1 + lp-ilp1)) * (q(iqp)**(il1+ilp1)) * gg * sumls(ix)
                    end do
                 end do
              end do

           end do
        end do

        do iq=1,nq
           do iqp=1,nq
              do ix=1,nx
                 gfac(a,ap,iq,iqp,ix) = sumq(iq,iqp,ix)
              end do
           end do
        end do

     end do
  end do

end subroutine permutation
!================================================================================================================================!
subroutine tsingle(ch,nptot,np,e,p2,wp,v,t)
use nrtype
use constants
implicit none
!--------------------------------------------------------------------!
! This subroutine determines the off-shell two-body t matrix
! at the energy e
! The subroutine dgesv from the lapack library solves the
! algebraic system:  akern * t = v
!--------------------------------------------------------------------!
integer(kind=i4b), intent(in) :: ch,nptot,np
real(kind=dp), intent(in) :: e
real(kind=dp), dimension(nptot), intent(in) :: p2,wp
real(kind=dp), dimension(6,nptot,nptot), intent(in) :: v
real(kind=dp), dimension(np,np), intent(out) :: t
integer(kind=i4b) :: i,ip,ifail
integer(kind=i4b), dimension(nptot) :: ipiv
real(kind=dp), dimension(nptot) :: wkspce
real(kind=dp), dimension(nptot,nptot) :: vv,akern

! Set up the potential

  do i=1,nptot
     do ip=1,nptot
        vv(i,ip)=v(ch,i,ip)
     end do
  end do

! Set up the t matrix kernel

  do i=1,nptot
     do ip=1,nptot
        akern(i,ip)=-vv(i,ip)*wp(ip)*p2(ip)/(e-(p2(ip)/mf))
        if(i == ip) akern(i,i)=akern(i,i)+one
     end do
  end do

! Solve the algebraic system

  ifail=4
  call dgesv(nptot,nptot,akern,nptot,ipiv,vv,nptot,ifail)
  if (ifail /= 0) then
     write(*,*) 'Error with dgsev in tmat'
     stop
  end if

  t(1:np,1:np)=vv(1:np,1:np)

end subroutine tsingle
!================================================================================================================================!
subroutine tcoupled(nptot,np,e,p2,wp,v,tpp,tmm,tpm,tmp)
use nrtype
use constants
implicit none
!--------------------------------------------------------------------!
! This subroutine determines the off-shell two-body t matrix
! at the energy e
! The subroutine dgesv from the lapack library solves the
! algebraic system:  akern * t = v
!--------------------------------------------------------------------!
integer(kind=i4b), intent(in) :: nptot,np
real(kind=dp), intent(in) :: e
real(kind=dp), dimension(nptot), intent(in) :: p2,wp
real(kind=dp), dimension(6,nptot,nptot), intent(in) :: v
real(kind=dp), dimension(np,np), intent(out) :: tpp,tmm,tpm,tmp
integer(kind=i4b) :: i,ip,ifail
integer(kind=i4b), dimension(2*nptot) :: ipiv
real(kind=dp), dimension(2*nptot) :: wkspce
real(kind=dp), dimension(nptot,nptot) :: vpp,vmm,vpm,vmp,identity,app,amm,apm,amp
real(kind=dp), dimension(2*nptot,2*nptot) :: vv,akern

! Set up the potential

  do i=1,nptot
     do ip=1,nptot
        vpp(i,ip) = v(3,i,ip)
        vmm(i,ip) = v(4,i,ip)
        vpm(i,ip) = v(5,i,ip)
        vmp(i,ip) = v(6,i,ip)
     end do
  end do

! Set up the identity matrix

  identity=zero
  do i=1,nptot
     identity(i,i)=one
  end do

  do i=1,nptot
     do ip=1,nptot
        app(i,ip)=-vpp(i,ip)*wp(ip)*p2(ip)/(e-p2(ip)/mf)
        amm(i,ip)=-vmm(i,ip)*wp(ip)*p2(ip)/(e-p2(ip)/mf)
        apm(i,ip)=-vpm(i,ip)*wp(ip)*p2(ip)/(e-p2(ip)/mf)
        amp(i,ip)=-vmp(i,ip)*wp(ip)*p2(ip)/(e-p2(ip)/mf)
     end do
  end do
  app=identity+app
  amm=identity+amm

! Set up the t matrix kernel and the interaction matrix

  vv(1:nptot,1:nptot)=vpp(1:nptot,1:nptot)
  vv(1:nptot,nptot+1:2*nptot)=vpm(1:nptot,1:nptot)
  vv(nptot+1:2*nptot,1:nptot)=vmp(1:nptot,1:nptot)
  vv(nptot+1:2*nptot,nptot+1:2*nptot)=vmm(1:nptot,1:nptot)

  akern(1:nptot,1:nptot)=app(1:nptot,1:nptot)
  akern(1:nptot,nptot+1:2*nptot)=apm(1:nptot,1:nptot)
  akern(nptot+1:2*nptot,1:nptot)=amp(1:nptot,1:nptot)
  akern(nptot+1:2*nptot,nptot+1:2*nptot)=amm(1:nptot,1:nptot)

! Solve the algebraic system

  ifail=4
  call dgesv(2*nptot,2*nptot,akern,2*nptot,ipiv,vv,2*nptot,ifail)
  if (ifail /= 0) then
     write(*,*) 'Error with dgsev in tmat'
     stop
  end if
  tpp(1:np,1:np)=vv(1:np,1:np)
  tpm(1:np,1:np)=vv(1:np,nptot+1:nptot+1+np)
  tmp(1:np,1:np)=vv(nptot+1:nptot+1+np,1:np)
  tmm(1:np,1:np)=vv(nptot+1:nptot+1+np,nptot+1:nptot+1+np)

end subroutine tcoupled
!================================================================================================================================!
subroutine search(f1,f2,e1,e2,e3,de)
use nrtype
use constants
implicit none
!---------------------------------------------------------------------------!
!  This subroutine is a simple search routine for a zero of a function
!---------------------------------------------------------------------------!
real(kind=dp), intent(in) :: f1,f2,e1,e2,de
real(kind=dp), intent(out) :: e3
integer(kind=i4b), save :: iv
real(kind=dp) :: v1,v2

  if (iv == 1) then
     e3=e2-f2*(e2-e1)/(f2-f1)
     iv=1
  else
     v1=sign(1.0_dp,f1)
     v2=sign(1.0_dp,f2)
     if (v1 /= v2) then
        e3=e2-f2*(e2-e1)/(f2-f1)
        iv=1
     else if (v1 == v2) then
        e3=e2+de/hbarc
        iv=0
     end if
  end if

end subroutine search
!================================================================================================================================!
end module trimod_sub