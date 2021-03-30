!*******************************************************
!*    PDMA algorithm subroutine for F90                *
!*                                                     *
!*       D. Hindle, University of Göttingen            *
!*                                                     *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* Simone Sebben & B. Rabi Baliga (1995)               * 
!* SOME EXTENSIONS OF TRIDIAGONAL AND PENTADIAGONAL    *
!* MATRIX ALGORITHMS, Numerical Heat Transfer, Part B: * 
!* Fundamentals: An International Journal of           *
!* Computation and Methodology,28:3, 323-351,          * 
!* DOI:10.1080/10407799508928837                       *
!*                                                     *
!*                                                     * 
!*******************************************************
MODULE PDMA

CONTAINS

!  ***************************************************************
! Given an N x N matrix A, which is pentadiagonal, 
! which is equivalent to a linear system
! au(i) + bu(i+1) + cu(i+2) + du(i-1) + eu(i-2) = f(i)
! where a = main diagonal, d and e are lr, b and c are ur diagonals and f is the load vector
! after rearranging for PDMA recursive algorithm 
! au(i) = f(i) - (bu(i+1) + cu(i+2) + du(i-1) + eu(i-2))  hence, carry out -1*(b+c+d+e)
! store matrix entries in vectors a(N), b(N), c(N), d(N), e(N), f(N)
! where d(1), e(1), e(2), b(N-1), c(N-1), c(N-2) lie outside the matrix Aij and are always zero
! and use recurrence relationship 
! 
! forwards from i=3,N (i=1,2 are first calculated directly)
! then set
! u(N) = r(N), u(N-1) = p(N-1) u(N) + r(N-1)  (these will be boundary conditions)
! implement backward recursion on
! u(i) = r(i) + p(i) u(i+1) + q(i) u(i+2) 
! from i=N-2,1
!  ***************************************************************

subroutine penta(n,a,b,c,d,e,f,u)

! a - main diagonal
! e,d - lr diagonals (e leftmost)
! b,c - ur diagonals (c rightmost)
! p,q,r, iteration variables
! f and u, vectors from system Au=f

! based on  Simone Sebben & B. Rabi Baliga (1995) SOME EXTENSIONS OF 
! TRIDIAGONAL AND PENTADIAGONAL MATRIX ALGORITHMS, Numerical Heat Transfer, Part B: 
! Fundamentals: An International Journal of Computation and Methodology,28:3, 323-351,
! DOI:10.1080/10407799508928837

! note, linear system is
! au(i) + bu(i+1) + cu(i+2) + du(i-1) + eu(i-2) = f(i)
! rearranging for PDMA recursive algorithm requires
! au(i) = f(i) - (bu(i+1) + cu(i+2) + du(i-1) + eu(i-2))  
! hence, to correspond to Sebben's algorithm, carry out -1*(b+c+d+e)

implicit none

integer, parameter :: QP = selected_real_kind (32)

real (kind=QP) :: p(n),q(n),r(n)
real (kind=QP) :: a(n),b(n),c(n),d(n),e(n),f(n),u(n)
real (kind=QP) :: bp(n),cp(n),dp(n),ep(n)
real (kind=QP) :: edpp(n),edpq(n),edpr(n),denom(n)
real (kind=QP), parameter :: zero=0.0_QP
integer :: i,j,k,l,n

!initialise
edpq=zero; edpp=zero; edpr=zero; denom=zero
p=zero; q=zero; r=zero; u=zero

! rearrange linear system as described in Sebben
do i = 1,n
  bp(i)=-b(i)
  cp(i)=-c(i)
  dp(i)=-d(i)
  ep(i)=-e(i)
end do

!set known iteration variable values for i=1, i=2, (p,q,r)

p(1)= bp(1) / a(1)  ! = 0
q(1)= cp(1) / a(1)  ! = 0
r(1)= f(1) / a(1)  ! = 0

p(2)= (bp(2) + dp(2)*q(1)) / (a(2) - dp(2)*p(1))  ! = 0
q(2)= cp(2) / (a(2) - dp(2)*p(1))                ! = 0  
r(2)= (f(2) + dp(2)*r(1))  / (a(2) - dp(2)*p(1)) ! = 0

! carry out forward recursion
! i=3,n 

do i=3,n
  edpq(i)=(dp(i)+ep(i)*p(i-2))*q(i-1)
  edpp(i)=(dp(i)+ep(i)*p(i-2))*p(i-1)
  edpr(i)=(dp(i)+ep(i)*p(i-2))*r(i-1)
  denom(i)=a(i)-ep(i)*q(i-2)-edpp(i)
  p(i)= (bp(i) + edpq(i)) / denom(i)
  q(i)= cp(i) / denom(i)
  r(i)= (f(i) + ep(i)*r(i-2) + edpr(i)) / denom(i)
end do

! set end values of u (n, n-1) to start backward recursion

u(n)=r(n)
u(n-1)=p(n-1)*u(n)+r(n-1)

! backward recursion algorithm to solve for u (n-2,n-1,...,1)

do i=n-2,1,-1
  u(i)=r(i) + p(i)*u(i+1) + q(i)*u(i+2)
end do

return

END subroutine penta

!---------------------------------------------------------------------------------------------
! subroutine to multiply a pentadiagonal matrix and a vector (only elements of rows i=3, node-2)
! u = a*f

subroutine pdmm(n,a1,am1,am2,ap1,ap2,f,u)

implicit none

integer, parameter :: QP = selected_real_kind (32)

real (kind=QP) :: a1(n),am1(n),am2(n),ap1(n),ap2(n),f(n),u(n)
integer :: i,j,k,l,n
u(1:n)=0.0_QP

do i=3,n-2
  u(i) = am2(i)*f(i-2)+am1(i)*f(i-1)+a1(i)*f(i)+ap1(i)*f(i+1)+ap2(i)*f(i+2)     
end do

return

END subroutine pdmm

END MODULE PDMA

! end of file pdma.f90

