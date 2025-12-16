MODULE nr_sub
IMPLICIT NONE
PRIVATE
PUBLIC :: gauleg,plgndr_s

  INTERFACE gauleg
          MODULE PROCEDURE gauleg
  END INTERFACE

  INTERFACE plgndr_s
          MODULE PROCEDURE plgndr_s
  END INTERFACE

CONTAINS
!================================================================================================================================!
SUBROUTINE gauleg(x1,x2,x,w)
USE nrtype; USE nrutil, ONLY : arth,assert_eq,nrerror
IMPLICIT NONE
REAL(DP), INTENT(IN) :: x1,x2
REAL(DP), DIMENSION(:), INTENT(OUT) :: x,w
REAL(DP), PARAMETER :: EPS=3.0e-14_dp
INTEGER(I4B) :: its,j,m,n
INTEGER(I4B), PARAMETER :: MAXIT=10
REAL(DP) :: xl,xm
REAL(DP), DIMENSION((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
LOGICAL(LGT), DIMENSION((size(x)+1)/2) :: unfinished
n=assert_eq(size(x),size(w),'gauleg')
m=(n+1)/2
xm=0.5_dp*(x2+x1)
xl=0.5_dp*(x2-x1)
z=cos(PI_D*(arth(1,1,m)-0.25_dp)/(n+0.5_dp))
unfinished=.true.
do its=1,MAXIT
	where (unfinished)
		p1=1.0
		p2=0.0
	end where
	do j=1,n
		where (unfinished)
			p3=p2
			p2=p1
			p1=((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
		end where
	end do
	where (unfinished)
		pp=n*(z*p1-p2)/(z*z-1.0_dp)
		z1=z
		z=z1-p1/pp
		unfinished=(abs(z-z1) > EPS)
	end where
	if (.not. any(unfinished)) exit
end do
if (its == MAXIT+1) call nrerror('too many iterations in gauleg')
x(1:m)=xm-xl*z
x(n:n-m+1:-1)=xm+xl*z
w(1:m)=2.0_dp*xl/((1.0_dp-z**2)*pp**2)
w(n:n-m+1:-1)=w(1:m)
END SUBROUTINE gauleg
!================================================================================================================================!
FUNCTION plgndr_s(l,m,x)
USE nrtype; USE nrutil, ONLY : arth,assert
IMPLICIT NONE
INTEGER(I4B), INTENT(IN) :: l,m
REAL(DP), INTENT(IN) :: x
REAL(DP), INTRINSIC :: abs,sqrt,product
REAL(DP) :: plgndr_s
INTEGER(I4B) :: ll
REAL(DP) :: pll,pmm,pmmp1,somx2

! Computes the associated Legendre polynomial P^l_m (x). Here m and l are integers satisfying
! 0 ≤ m ≤ l, while x lies in the range −1 ≤ x ≤ 1.

  call assert(m >= 0, m <= l, abs(x) <= 1.0_dp, 'plgndr_s args')
  pmm=1.0_dp
  if (m > 0) then
     somx2=sqrt((1.0_dp-x)*(1.0_dp+x))
     pmm=product(arth(1.0_dp,2.0_dp,m))*somx2**m
     if (mod(m,2) == 1) pmm=-pmm
  end if
  if (l == m) then
     plgndr_s=pmm
  else
     pmmp1=x*(2*m+1)*pmm
     if (l == m+1) then
        plgndr_s=pmmp1
     else
        do ll=m+2,l
           pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
           pmm=pmmp1
           pmmp1=pll
        end do
        plgndr_s=pll
     end if
  end if

END FUNCTION plgndr_s
!================================================================================================================================!
END MODULE nr_sub
