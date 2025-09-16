module trimod_sub
implicit none
private
public :: gauss,pot,tmat,search

  interface gauss
          module procedure gauss
  end interface

  interface pot
          module procedure pot
  end interface

  interface tmat
          module procedure tmat
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

  if (n >= 50) then
     write(*,*) 'Error in gauss: n > 49'
     stop
  end if

  call gauleg(-1.0_dp,1.0_dp,gp,gw)
  b1=2.*a*c-a*b-b*c
  b2=2.*a*(b-c)
  b3=a-2.*b+c
  b4=a-c
  b5=b-c
  b6=b-a
  do i=1,n
     cn=b3*gp(i)+b4
     x(i)=(b1*(1.+gp(i))+b2)/cn
     x2(i)=x(i)*x(i)
     wx(i)=x2(i)*2.*b5*b4*b6*gw(i)/cn**2
  end do

end subroutine gauss
!================================================================================================================================!
subroutine pot(nptot,p,p2,wp,v,v0,a)
use nrtype
use constants
implicit none
!--------------------------------------------------------------------!
! This subroutine determines a Malfliet-Tjion potential
! in momentum space
!--------------------------------------------------------------------!
integer(kind=i4b), intent(in) :: nptot
real(kind=dp), dimension(nptot), intent(in) :: p,p2,wp
real(kind=dp), dimension(nptot,nptot), intent(out) :: v
integer(kind=i4b) :: i,ip,ipp
real(kind=dp) :: x1,x2,vv,qv,z,x
!real(kind=dp), dimension(2) :: v0,a,v0m
! AGGIUNGO I POTENZIALI PASSATI DA ESTERNO:
real(kind=dp), dimension(2), intent(in) :: v0, a
real(kind=dp), dimension(2) :: v0m

!  v0(1)=-570.316_dp
!  v0(2)=1438.4812_dp
!  a(1)=1.55_dp
!  a(2)=3.11_dp

  do i=1,2
     v0m(i)=(v0(i)/hbarc)*mf
  end do

! Set up potential

  do ip=1,nptot
     do ipp=1,ip
        x1=p(ip)
        x2=p(ipp)
        vv=0.0_dp
        do i=1,2
           if (x1*x2 == 0.0_dp) then
              qv=2.0_dp/(x1**2+x2**2+a(i)**2)
           else
              z=(a(i)**2+x1**2+x2**2)/(2.0_dp*x1*x2)
              if (z < 1.0e2) then
                 qv=(log(z+1.0_dp) - log(z-1.0_dp))/(2.0_dp*x1*x2)
              else
                 x=1.0_dp/z
                 qv=(x+x**3/3.0_dp+x**5/5.0_dp+x**7/7.0_dp)*2.0_dp/(2.0_dp*x1*x2)
              end if
           end if
           vv=vv+qv*v0m(i)/pi_d
        end do
        v(ip,ipp)=vv
        v(ipp,ip)=v(ip,ipp)
     end do
  end do

end subroutine pot
!================================================================================================================================!
subroutine tmat(e,nptot,p,p2,wp,v,t)
use nrtype
implicit none
!--------------------------------------------------------------------!
! This subroutine determines the off-shell two-body t matrix
! at the energy e
! The subroutine dgesv from the lapack library solves the
! algebraic system:  akern * t = v
!--------------------------------------------------------------------!
integer(kind=i4b), intent(in) :: nptot
real(kind=dp), intent(in) :: e
real(kind=dp), dimension(nptot), intent(in) :: p,p2,wp
real(kind=dp), dimension(nptot,nptot), intent(in) :: v
real(kind=dp), dimension(nptot,nptot), intent(out) :: t
integer(kind=i4b) :: i,ip,ifail
integer(kind=i4b), dimension(nptot) :: ipiv
real(kind=dp), dimension(nptot) :: wkspce
real(kind=dp), dimension(nptot,nptot) :: akern,vv

! Set up the t matrix kernel

  do i=1,nptot
     do ip=1,nptot
        akern(i,ip)=-v(i,ip)*wp(ip)/(e-p2(ip))
        if (i == ip) akern(i,i)=akern(i,i)+1.0_dp
     end do
  end do

! Solve the algebraic system

  ifail=4
  vv=v
  call dgesv(nptot,nptot,akern,nptot,ipiv,vv,nptot,ifail)
  if (ifail /= 0) then
     write(*,*) 'Error with dgsev in tmat'
     stop
  end if
  t=vv

end subroutine tmat
!================================================================================================================================!
subroutine search(f1,f2,e1,e2,e3)
use nrtype
use constants
implicit none
!---------------------------------------------------------------------------!
!  This subroutine is a simple search routine for a zero of a function
!---------------------------------------------------------------------------!
real(kind=dp), intent(in) :: f1,f2,e1,e2
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
        e3=e2+(de/hbarc)*mf
        iv=0
     end if
  end if

end subroutine search
!================================================================================================================================!
end module trimod_sub