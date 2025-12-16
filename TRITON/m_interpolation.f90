module interpolation
implicit none
private
public :: spline

  interface spline
          module procedure spline
  end interface

contains
!================================================================================================================================!
subroutine spline(nptot,np,nq,nx,q,q2,p,x,wx,s1,s2)
use nrtype
use nr_sub
implicit none
!--------------------------------------------------------------------!
! This subroutine determines the spline elements
!--------------------------------------------------------------------!
integer(kind=i4b), intent(in) :: nptot,np,nq,nx
real(kind=dp), dimension(nq), intent(in) :: q,q2
real(kind=dp), dimension(nptot), intent(in) :: p
real(kind=dp), dimension(nx), intent(out) :: x,wx
real(kind=dp), dimension(nx,np,nq,nq), intent(out) :: s1,s2
integer(kind=i4b) :: iq,iqp,ix,ip
real(kind=dp) :: c1,c2,c3,c4,pi1,pi2
real(kind=dp), dimension(np) :: spl1,spl2
real(kind=dp), dimension(40,40) :: fak1,fak2,fak3

  call sprep(np,p,fak1,fak2,fak3)

  call gauleg(-1.0_dp,1.0_dp,x,wx)

  do iq=1,nq
     do iqp=1,nq
        c1=0.25_dp*q2(iq)+q2(iqp)
        c2=q2(iq)+0.25_dp*q2(iqp)
        c3=q(iq)*q(iqp)
        do ix=1,nx
           c4=c3*x(ix)
           pi1=sqrt(c1+c4)
           pi2=sqrt(c2+c4)
           call selem(np,p,fak1,fak2,fak3,pi1,spl1)
           call selem(np,p,fak1,fak2,fak3,pi2,spl2)
           do ip=1,np
              s1(ix,ip,iqp,iq)=spl1(ip)
              s2(ix,ip,iqp,iq)=spl2(ip)
           end do
        end do
     end do
  end do

end subroutine spline
!================================================================================================================================!
subroutine sprep(n,x,fak1,fak2,fak3)
use nrtype
implicit none
!--------------------------------------------------------------------!
! This subroutine prepares coefficients for the spline elements
!--------------------------------------------------------------------!
integer(kind=i4b), intent(in) :: n
real(kind=dp), dimension(n), intent(in) :: x
real(kind=dp), dimension(40,40), intent(out) :: fak1,fak2,fak3
integer(kind=i4b) :: i,j
real(kind=dp) :: ax,bx,cx,al,am,pim,h1
real(kind=dp), dimension(40) :: hi,u
real(kind=dp), dimension(40,40) :: q,c

  if (n > 40) stop
  u(1)=0.0_dp
  hi(2)=x(2)-x(1)
  do i=1,n
     q(1,i)=0.0_dp
  end do
  do i=2,n-1
     ax=x(i+1)-x(i)
     hi(i+1)=ax
     bx=x(i+1)-x(i-1)
     cx=x(i)-x(i-1)
     al=ax/bx
     am=1.0_dp-al
     pim=1.0_dp/(2.0_dp-am*u(i-1))
     u(i)=al*pim
     do j=1,n
        q(i,j)=-pim*am*q(i-1,j)
     end do
     q(i,i-1)=q(i,i-1)+pim/(cx*bx)
     q(i,i)=q(i,i)-pim/(cx*ax)
     q(i,i+1)=q(i,i+1)+pim/(ax*bx)
  end do
  do j=1,n
     c(n,j)=0.0_dp
     fak1(N,J)=0.0_dp
     fak2(N,J)=0.0_dp
     fak3(N,J)=0.0_dp
  end do
  do i=n-1,1,-1
     h1=1.0_dp/hi(i+1)
     do j=1,n
        c(i,j)=q(i,j)-c(i+1,j)*u(i)
        fak1(i,j)=-hi(i+1)*(2.0_dp*c(i,j)+c(i+1,j))
     end do
     fak1(i,i)=fak1(i,i)-h1
     fak1(i,i+1)=fak1(i,i+1)+h1
     do j=1,n
        fak2(i,j)=3*c(i,j)
        fak3(i,j)=(c(i+1,j)-c(i,j))*h1
     end do
  end do

end subroutine sprep
!================================================================================================================================!
subroutine selem(n,x,fak1,fak2,fak3,xa,spl)
use nrtype
implicit none
!--------------------------------------------------------------------!
! This subroutine determines the spline elements spl
!--------------------------------------------------------------------!
integer(kind=i4b), intent(in) :: n
real(kind=dp), dimension(n), intent(in) :: x
real(kind=dp), dimension(40,40), intent(in) :: fak1,fak2,fak3
real(kind=dp), intent(in) :: xa
real(kind=dp), dimension(n), intent(out) :: spl
integer(kind=i4b) :: i,j
real(kind=dp) :: dx

  if (xa < x(1) .or. xa > x(n)) then
     write(*,*) 'Error in selem'
     stop
  end if
  i=-1
  do
     i=i+1
     if (xa >= x(i+1) .and. i <= n) then
        cycle
     else
        exit
     end if
  end do
  if (i == 0) i=1
  dx=xa-x(i)
  do j=1,n
     spl(j)=((fak3(i,j)*dx+fak2(i,j))*dx+fak1(i,j))*dx
  end do
  spl(i)=spl(i)+1.0_dp

end subroutine selem
!================================================================================================================================!
end module interpolation
