module trimod_sub
implicit none
private
public :: gauss,pot,tmat,search, tmat_chiral_ekstrom

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

  
interface tmat_chiral_ekstrom
    module procedure tmat_chiral_ekstrom
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
  use m_params,     only: ipot
  use m_pot_chiral, only: vchi, chiral_loaded, nptot_chi
  implicit none
  ! ------- dichiarazioni -------
  integer(kind=i4b), intent(in) :: nptot
  real(kind=dp), dimension(nptot), intent(in)  :: p,p2,wp
  real(kind=dp), dimension(nptot,nptot), intent(out) :: v
  real(kind=dp), dimension(2), intent(in) :: v0, a

  integer(kind=i4b) :: i,ip,ipp
  real(kind=dp) :: x1,x2,vv,qv,z,x
  real(kind=dp), dimension(2) :: v0m

  logical, save :: first = .true.
  integer :: upot, ios_p
  real(kind=dp), parameter :: alpha_chiral = 9.2e-9_dp
  ! -----------------------------------------------

    ! === Caso 2: potenziale chirale tabulato ===
  if (ipot == 2) then
     

     if (.not. chiral_loaded) then
        write(*,*) 'Error in pot: chiral potential not loaded!'
        stop
     end if
     if (nptot > nptot_chi) then
        write(*,*) 'Error in pot: nptot > nptot_chi'
        stop
     end if

     ! copia la matrice chiral SCALATA in v
     v(1:nptot,1:nptot) = alpha_chiral * vchi(1:nptot,1:nptot)

     write(*,'(A,2(1X,ES12.5))') 'DEBUG V (scaled) min,max =',  &
          minval(v(1:nptot,1:nptot)), maxval(v(1:nptot,1:nptot))

     return
  end if


  ! === Caso default: Malfliet–Tjon (potenziale originale) ===

  ! Fattori di conversione per il potenziale MT
  do i=1,2
     v0m(i) = (v0(i)/hbarc)*mf
  end do

  ! Blocco di debug: scrivi i parametri usati, una sola volta
  if (first) then
     open(newunit=upot, file='pot.used', status='replace', &
          action='write', iostat=ios_p)
     if (ios_p == 0) then
        write(upot,'(A,2(1X,ES16.8),A,2(1X,ES16.8))') &
             'v0=', v0(1), v0(2), 'a=', a(1), a(2)
        close(upot)
     end if
     first = .false.
  end if

  ! Costruzione del potenziale in spazio dei momenti (Malfliet–Tjon)
  do ip=1,nptot
     do ipp=1,ip
        x1 = p(ip)
        x2 = p(ipp)
        vv = 0.0_dp
        do i=1,2
           if (x1*x2 == 0.0_dp) then
              qv = 2.0_dp/(x1**2 + x2**2 + a(i)**2)
           else
              z = (a(i)**2 + x1**2 + x2**2)/(2.0_dp*x1*x2)
              if (z < 1.0e2_dp) then
                 qv = (log(z+1.0_dp) - log(z-1.0_dp))/(2.0_dp*x1*x2)
              else
                 x  = 1.0_dp/z
                 qv = (x + x**3/3.0_dp + x**5/5.0_dp + x**7/7.0_dp) * &
                       2.0_dp/(2.0_dp*x1*x2)
              end if
           end if
           vv = vv + qv*v0m(i)/pi_d
        end do
        v(ip,ipp) = vv
        v(ipp,ip) = vv
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





! nuova Ls "compatibile" con Ekstrom
!================================================================================================================================!
subroutine tmat_chiral_ekstrom(e2, nptot, p, v, t)
  use nrtype
  use constants, only: hbarc, m, pi_d
  implicit none
  !----------------------------------------------------------------
  ! LS 2-body in stile Ekström:
  !
  ! T = ( I - (2/pi) * V * G0(E) * 2*mu )^{-1} * V
  !
  ! dove:
  !   p       : momenta TRIMOD in fm^-1
  !   e2      : energia 2-body in unità interne TRIMOD (fm^-2)
  !   V       : matrice del potenziale chiral in unità "di Ekström"
  !             già letta da v_chiral_matrix.dat
  !   T       : matrice t off–shell in "unità LS Ekström"
  !
  ! Per now:
  !   p_MeV  = p * hbarc
  !   E2_MeV = e2 * hbarc^2 / m
  !   mu     = m/2
  !
  !----------------------------------------------------------------
  integer(kind=i4b), intent(in) :: nptot
  real(kind=dp),     intent(in) :: e2
  real(kind=dp),     dimension(nptot), intent(in)  :: p
  real(kind=dp),     dimension(nptot,nptot), intent(in)  :: v
  real(kind=dp),     dimension(nptot,nptot), intent(out) :: t

  ! workspace
  real(kind=dp), dimension(nptot)       :: p_MeV, wq, eps_k, G0
  real(kind=dp), dimension(nptot,nptot) :: K, A, VV
  integer(kind=i4b), dimension(nptot)   :: ipiv
  integer(kind=i4b) :: i,j
  real(kind=dp) :: mu, E2_MeV, normw

  !------------------------
  ! 1) converto p in MeV
  !------------------------
  do i = 1, nptot
     p_MeV(i) = p(i) * hbarc      ! [MeV]
  end do

  ! massa ridotta NN (approssimiamo con m/2)
  mu = 0.5_dp * m                 ! [MeV]

  ! energia relativa 2-body in MeV:
  ! e2 (fm^-2) -> E2_MeV = e2 * (hbarc^2/m)
  E2_MeV = e2 * (hbarc*hbarc) / m

  !------------------------
  ! 2) pesi di integrazione dq (trapezoidale su p_MeV)
  !    wq(i) ~ Δq_i in MeV
  !------------------------
  if (nptot >= 2) then
     wq(1) = 0.5_dp * (p_MeV(2) - p_MeV(1))
     do i = 2, nptot-1
        wq(i) = 0.5_dp * (p_MeV(i+1) - p_MeV(i-1))
     end do
     wq(nptot) = 0.5_dp * (p_MeV(nptot) - p_MeV(nptot-1))
  else
     wq(1) = 0.0_dp
  end if

  ! (facoltativo) normalizza un minimo per debug
  normw = 0.0_dp
  do i = 1, nptot
     normw = normw + wq(i)
  end do
  ! write(*,'(A,ES12.5)') 'DEBUG sum dq =', normw

  !------------------------
  ! 3) propagatore libero G0(E; p)
  !    G0(j) = 1 / (E2_MeV - p^2/(2 mu))
  !------------------------
  do j = 1, nptot
     eps_k(j) = (p_MeV(j)*p_MeV(j)) / (2.0_dp*mu)
     G0(j)    = E2_MeV - eps_k(j)
     if (abs(G0(j)) < 1.0e-6_dp) then
        G0(j) = sign(1.0e-6_dp, G0(j))
     end if
     G0(j) = 1.0_dp / G0(j)
  end do

  !------------------------
  ! 4) Kernel K = (2/pi) * V * G0 * (2 mu) * q^2 dq
  !    discretizzazione:
  !    K(i,j) = (2/pi) * V(i,j) * (2*mu) * G0(j) * p_MeV(j)^2 * wq(j)
  !------------------------
  do i = 1, nptot
     do j = 1, nptot
        K(i,j) = (2.0_dp/pi_d) * v(i,j) * (2.0_dp*mu) * G0(j) &
                 * (p_MeV(j)*p_MeV(j)) * wq(j)
     end do
  end do

  !------------------------
  ! 5) Risolvi (I - K) T = V
  !------------------------
  do i = 1, nptot
     do j = 1, nptot
        A(i,j) = -K(i,j)
     end do
     A(i,i) = A(i,i) + 1.0_dp
  end do

  VV = v
  call dgesv(nptot, nptot, A, nptot, ipiv, VV, nptot, i)

  if (i /= 0) then
     write(*,*) 'Error in dgesv inside tmat_chiral_ekstrom, info =', i
     stop
  end if

  t = VV

end subroutine tmat_chiral_ekstrom





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
integer(kind=i4b), save :: iv = 0
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