module mod_factorial
implicit none
integer, parameter :: ndfact=70
integer :: nf
real(8) :: xafac(0:ndfact)

end module mod_factorial
!================================================================================================================================!
subroutine factorial
use mod_factorial
implicit none
integer :: i

! Factorials for xg, threej,... functions

  nf=70
  xafac(0)=1.0d0
  do i=1,nf
     xafac(i)=real(i)*xafac(i-1)*0.1d0
  end do

  return
end subroutine factorial
!================================================================================================================================!
function fact(n)
implicit none
real(8) :: fact
integer :: n, i
integer(8) :: nfact

! Calculate n!

  nfact=1
  do i=1,n
     nfact = nfact*i
  end do
  
  fact=dfloat(nfact)
 
  return
end function fact
!================================================================================================================================!
real(8) function xg(l1,m1,l2,m2,l,m)
implicit real(8) (a-h,o-z)
!----------------------------------------
! Clebsch-Gordan coefficients
! N.B. l1,m1,l2,m2,l,m are doubled
!----------------------------------------
  xg=0.d0
  if (m1+m2 .ne. m) return
  if (l .lt. iabs(l1-l2) .or. l .gt. (l1+l2)) return
  if (iabs(m1) .gt. l1) return
  if (iabs(m2) .gt. l2) return
  if (iabs(m) .gt. l) return

  ll=(l1-l2+m)/2
  sign=real((-1)**ll)
  xg=sign*sqrt(real(l+1))*threej(l1,l2,l,m1,m2,-m)
  
  return
end function xg
!================================================================================================================================!
real(8) function threej(j1,j2,j,m1,m2,m)
use mod_factorial
!----------------------------------------
! Three-j symbol
! N.B. j's and m's ====> 2j's and 2m's
!      (j1+j2+j)/2+1 < nf
!----------------------------------------
implicit real(8) (a-h,o-z)
integer :: z,zmin,zmax

    cc=0.0d0
    if (m1+m2 .ne. -m .or. iabs(m1) .gt. iabs(j1) .or. iabs(m2) .gt. iabs(j2) &
       .or. iabs(m) .gt. iabs(j) .or. j .gt. j1+j2 .or. j .lt. iabs(j1-j2)) goto 1
    zmin=0
    if (j-j2+m1 .lt. 0) zmin=-j+j2-m1
    if (j-j1-m2+zmin .lt. 0) zmin=-j+j1+m2
    zmax=j1+j2-j
    if (j2+m2-zmax .lt. 0) zmax=j2+m2
    if (j1-m1-zmax .lt. 0) zmax=j1-m1
    do z=zmin,zmax,2
       ja=z/2
       jb=(j1+j2-j-z)/2
       jc=(j1-m1-z)/2
       jd=(j2+m2-z)/2
       je=(j-j2+m1+z)/2
       jf=(j-j1-m2+z)/2
       fase=real((-1)**(z/2))
       cc=cc+fase/(xafac(ja)*xafac(jb)*xafac(jc)*xafac(jd)*xafac(je)*xafac(jf))
    end do
    ja=(j1+j2-j)/2
    jb=(j1-j2+j)/2
    jc=(-j1+j2+j)/2
    jd=(j1+m1)/2
    je=(j1-m1)/2
    jf=(j2+m2)/2
    jg=(j2-m2)/2
    jh=(j-m)/2
    ji=(j+m)/2
    jj=(j1+j2+j+2)/2
    if (jj .gt. nf) call quit(2,'Error - Factorials exceeded')
    cc=sqrt(real(j+1)*xafac(ja)*xafac(jb)*xafac(jc)*xafac(jd))*cc
    cc=sqrt(xafac(je)*xafac(jf)*xafac(jg)*xafac(jh)*xafac(ji)/xafac(jj))*cc

1   threej=real((-1)**((j1-j2-m)/2))*cc/sqrt(10.0d0*real(j+1))

    return
end function threej
!================================================================================================================================!
real(8) function sixj(j1,j2,j3,l1,l2,l3)
use mod_factorial
!----------------------------------------
! Six-j symbol
! N.B. j's and l's ====> 2j's and 2l's
!----------------------------------------
implicit real(8) (a-h,o-z)
integer z,zmin,zmax

    cc=0.0d0
    d1=dtr(j1,j2,j3)
    d2=dtr(j1,l2,l3)
    d3=dtr(l1,j2,l3)
    d4=dtr(l1,l2,j3)
    if (d1 .eq. 0.0d0 .or. d2 .eq. 0.0d0 .or. d3 .eq. 0.0d0 .or. d4 .eq. 0.0d0) goto 1
    zmin=j1+j2+j3
    if (zmin-j1-l2-l3 .lt. 0) zmin=j1+l2+l3
    if (zmin-l1-j2-l3 .lt. 0) zmin=l1+j2+l3
    if (zmin-l1-l2-j3 .lt. 0) zmin=l1+l2+j3
    zmax=j1+j2+l1+l2
    if (j2+j3+l2+l3-zmax .lt. 0) zmax=j2+j3+l2+l3
    if (j3+j1+l3+l1-zmax .lt. 0) zmax=j3+j1+l3+l1
    do z=zmin,zmax,2
       ja=z/2+1
       jb=(z-j1-j2-j3)/2
       jc=(z-j1-l2-l3)/2
       jd=(z-l1-j2-l3)/2
       je=(z-l1-l2-j3)/2
       jf=(j1+j2+l1+l2-z)/2
       jg=(j2+j3+l2+l3-z)/2
       jh=(j3+j1+l3+l1-z)/2
       fase=real((-1)**(z/2))
       ccd=fase*xafac(ja)/xafac(jb)/xafac(jc)/xafac(jd)/xafac(je)/xafac(jf)
       cc=cc+ccd/xafac(jg)/xafac(jh)
    end do
1   sixj=d1*d2*d3*d4*cc*10.0d0

    return
end function sixj
!================================================================================================================================!
real(8) function dtr(j1,j2,j3)
use mod_factorial
!----------------------------------------
!  (j1+j2+j)/2+1 < nf
!----------------------------------------
implicit real(8) (a-h,o-z)

    if (iabs(j1-j2) .gt. j3 .or. j1+j2 .lt. j3) then
       dtr=0.0d0
    else
       jd=(j1+j2+j3)/2+1
       if (jd .gt. nf) call quit(2,'ERROR - FACTORIALS EXCEEDED')
       ja=(j1+j2-j3)/2
       jb=(j1-j2+j3)/2
       jc=(-j1+j2+j3)/2
       dtr=sqrt(xafac(ja)*xafac(jb)*xafac(jc)/xafac(jd)/10.0d0)
    end if

    return
end function dtr
!================================================================================================================================!
real(8) function aninej(j1,j2,j12,j3,j4,j34,j13,j24,j)
use mod_factorial
!----------------------------------------
! Nine-j symbol
! N.B. j's and l's ====> 2j's and 2l's
!----------------------------------------
implicit real(8) (a-h,o-z)

    cc=0.d0
    do jp=iabs(j-j1),j+j1,2
       cg1=sixj(j1,j3,j13,j24,j,jp)
       cg2=sixj(j2,j4,j24,j3,jp,j34)
       cg3=sixj(j12,j34,j,jp,j1,j2)
       si=real((-1)**jp)
       cc=cc+si*real(jp+1)*cg1*cg2*cg3
    end do
    aninej=cc

    return
end function aninej
!================================================================================================================================!
logical function triang(j1,j2,j3)
implicit none
integer j1,j2,j3

    triang=.false.
    if (j1 .gt. (j2+j3)) return
    if (j1 .lt. abs(j2-j3)) return
    if (mod(j1+j2+j3,2) .ne. 0) return
    triang=.true.

end function triang
!================================================================================================================================!
subroutine quit(mode,notice)
implicit none
!===========================================C
!  RUN TERMINATION WITH ERROR MESSAGE      C
!===========================================C
integer, parameter :: masof=6
character(17) :: remark(2)
character     :: notice*(*)
integer       :: mode
data remark /'INPUT ERROR','RUN TIME ERROR'/

    write (masof,'(//1x,a/1x,a/1x,a,a//)') &
          'EXECUTION IS TERMINATED BECAUSE AN ERROR IS FOUND>', &
          remark(MODE), 'DETAILS> ', NOTICE
    stop 'QUIT - ERROR'

end subroutine quit