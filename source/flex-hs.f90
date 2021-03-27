! program flexdqup
! calculate flexure of plate under dynamic sediment load - solved iteratively - - pinned at both ends
! solve 
! (D u'')'' + Pu'' + k u = q(u)  ; u(0)=a, u'(0)= b, u'(L)=c, u(L)=d
! boundary conditions fix plate at both ends
!--------------------------------------------------------------------------------------------
! PARAMETERS
! D = flexural rigidity = Eh^3/(12(1-v^2)) 
! EM = elastic modulus of elastic lithosphere (~7e10 Pa)
! h = effective elastic thickness of plate (~3e4 m)
! P = horizontal ("buckling") force from the ends of the plate (Nm, +ve = compression, reasonable value ~10^12 - 10^13)
! v = poisson's ratio (0.25)
! k=(pm)g (pm=3300 kgm-3, pw=1000 kgm-3)
! ps = load density (e.g. 2300 kgm-3)
! dx = grid spacing for numerical solution (5000 metres)
! -------------------------------------------------------------------------------------------
! METHOD
! 
! equation reduces to a finite difference approximation.
! D_i+1/h4 u_i+2 - 2(D_i+1 + D_i)/h4 u_i+1 + ((D_i+1 +4D_i + D_i-1)/h4 + k) u_i - 2(D_i-1 + D_i)/h4 u_i-1 + D_i-1/h4 u_i-2 = q_i
! leads to linear system
! a u = q(u)
! solve using pdma recursive algorithm with forward then back substitution
! needs an iterative solution for the changing value of q(u) 
! solve A-1 qu(r) = u(r+1) 
! solve A-1 qu(r+1) = u(r+2) etc. etc. until convergence
! Where u(r) u(r+1) refers to iteration step. 
!
!-------------------------------------------------------------------------------------------
! INPUT VARIABLES
! q(i) - sediment thickness at node i of a plate (lp, i=100)
! te(i) - elastic thickness at node i (wherever it deviates from zero)
! OUTPUT VARIABLES 	 
! u(i) - plate geometry at node i of a plate (lp, i=100)
!
!_____________________________________________________________________________________________
! ADDITIONAL VARIABLES NEEDED
! 1. convergence parameter/tolerance
! 2. uo, ur, ur1: first, preceding, and current solution to problem respectively
! ur1(i) - ur(i)  gives convergence parameter (sum_i ABS(u_i(r+1)-u_i(r)))
! 3. load modification code as a function of u(r+1) to calculate q(u(r+1))
! 4. iteration loop to make repeated calculations as necessary
!
! load subroutine module - PDMA for recursive pentadiagonal matrix algorithm for solving A^-1 q = u
! load subroutine module - PLOTFLEX for plot output (works with GMT v6 script plotflex6.gmt
 use PDMA
 use PLOTFLEX
! ------------------------------------------------

! declare variables

implicit none
integer, parameter :: QP = selected_real_kind (32)

integer :: node       ! size/resolution of the problem i.e. how many nodes across the domain, last element of line 3 in param.txt  
real (kind=QP), dimension (:), allocatable :: a1,ap1,ap2,am1,am2,q,u,t,tee,ur,ur1,uo,lf,vd,pu  
real (kind=QP) :: EM, h, v, k, dx, pm, pfill, pload, pcrust, g, pup, tol ! normal dp variables from param.txt
real (kind=QP) ::  te, plotf, th, udiff
real (kind=QP) :: bc1, bc2, bc3, bc4  ! 4 boundary condition variables
real (kind=QP) :: depth, pi,dum1,dum2,dum3
double precision, parameter :: zero=0.0_QP
integer :: i,j,l,m,n,r,iter,check,istore,il,ir,pint
character*3 :: code
logical :: iswitch,bc3rd,comp

code='hs3'
! open input AND output files
! INPUT   i.e these files must already exist and be correctly formatted - i.e contain the correct number and type of variables in the right places  
open (10, file='te.txt')        ! elastic thickness profile.... (h is still the default elastic thickness, read from param.txt)
open (11, file='load-flex.txt') ! load input for problem 
open (12, file='bcond.txt')     ! boundary conditions file
open (13, file='param.txt')     ! model parameters file
! OUTPUT

! Begin reading input....
! read parameters from param.txt 

open (13, file='param.txt')     ! model parameters file
read (13,*)
read (13,*)
read (13,*) EM, h, pm, pload, pfill, pcrust, g, dx, v, pup, tol, node, iswitch, pint, bc3rd, comp

! dynamic allocation of arrays, size given as parameter node in param.txt

allocate (a1(node),ap1(node),ap2(node),am1(node),am2(node),q(node),pu(node))

allocate (vd(node),tee(node),t(node),u(node),ur(node),uo(node),ur1(node),lf(node))

print *, 'allocated ', node, ' spaces in arrays'

k=g*pm

! initialise vectors
a1=zero; ap1=zero; ap2=zero; am1=zero; am2=zero; q=zero
vd=zero; t=zero; lf=zero; ur=zero; ur1=zero; uo=zero; pu=pup

tee=h

print *, 'initialised arrays'

! read in elastic thickness at node (i) (in metres)
l = 0
30 continue
  l = l+1
  read (10,*, end=70) il, ir,te
  if (il < 1 .or. ir > node) then
    write(*,'(a,i6,a)') 'error in file "te.txt", line',l,', inconsistent indices'
    stop
  end if
  do i=il,ir
    tee(i)=te
    pu(i)=te/h*pup
    write(*,'(1pe15.5)') pu(i)
  end do      
goto 30

70 continue

! calculate the value of DD(i) - up to this point, tee(i) carries correct value for elastic thickness
vd(1:node)=(EM*tee(1:node)**3)/(12.0_QP*(1.0_QP-v**2))

!  read in load thickness at i (in metres)
l = 0
33 continue
  l = l+1
  read (11,*, end=99) il,ir,th
  if (il < 1 .or. ir > node) then
    write(*,'(a,i6,a)') 'error in file "load-flex.txt", line',l,', inconsistent indices'
    stop
  end if
  do i=il,ir
    t(i)=th 
    q(i)=t(i)*g*pload
  end do
goto 33

! read boundary conditions from bcond.txt
99  continue

if (bc3rd .eqv. .true.) then 
  read (12,*)
  read (12,*)
end if
read (12,*)     
read (12,*) bc1, bc2, bc3, bc4
q(1)=bc1
q(2)=bc2
q(node-1)=bc3
q(node)=bc4

if (bc3rd .eqv. .true.) then
  q(2)=0.5_QP*bc2/vd(2)
  print *, 'q(2) = ', q(2) 
end if           

100  continue 

! 5 diagonals of matrix
do i=3,node-2
  a1(i) = (1._QP/dx**4)*(vd(i+1)+4._QP*vd(i)+vd(i-1)) -k - 2._QP*pu(i)/dx**2   ! i
  !a1(i) = (1._QP/dx**4)*(vd(i+1)+4._QP*vd(i)+vd(i-1)) -k    ! i
  ap1(i) = (-2._QP/dx**4.)*(vd(i+1)+vd(i)) + pu(i)/dx**2                   ! i+1  ur diag
  ap2(i) = (1._QP/dx**4)*vd(i+1)                                         ! i+2    ur diag
  am1(i) = (-2._QP/dx**4)*(vd(i)+vd(i-1)) + pu(i)/dx**2                    ! i-1  lr diag
  am2(i) = (1._QP/dx**4)*vd(i-1)                                         ! i-2    lr diag
end do

!   boundary conditions
if (bc3rd .eqv. .true.) then
  print *, '3rd order bc '

!   u''(0)

  a1(1)=1._QP/dx**2
  ap1(1)=-2._QP/dx**2
  ap2(1)=1._QP/dx**2

!   u'''(0)

  am1(2)=-1._QP/dx**3
  a1(2)=3._QP/dx**3
  ap1(2)=-3._QP/dx**3
  ap2(2)=1._QP/dx**3

else
 print *, 'normal bc '

!   u(0)  

 a1(1)=1._QP

!   u'(0)   forward diff, 2nd order

  am1(2)=-1.5_QP/(1._QP*dx)
  a1(2)=2._QP/(1._QP*dx)
  ap1(2)=-0.5_QP/(1._QP*dx)  
!  am1(2)=-1._QP/dx    ! ur diag
!  a1(2)=1._QP/dx

end if

!   u'(L)   backward diff, 2nd order

ap1(node-1)=1.5_QP/(1._QP*dx)
a1(node-1)=-2._QP/(1._QP*dx)
am1(node-1)=0.5_QP/(1._QP*dx)       
!a1(node-1)=-1._QP/dx
!ap1(node-1)=1._QP/dx   ! lr diag

!   u(L) 

a1(node)=1._QP

print *, am1(2),a1(2),ap1(2) 
print *, ap1(node-1),a1(node-1),am1(node-1)
print *, q(2),q(node-1)

if (iswitch .eqv. .false.) then 
call penta(node,a1,ap1,ap2,am1,am2,q,ur1)      ! solution gets stored in ur1 - should always be u(r+1)
goto 199
end if       

! intitialise udiff (iteration convergence)
udiff=tol+1._QP
iter=1

do while (udiff>tol)

  udiff=0._QP

! calculate A-1 q_r = u_r+1 - and iterate

  call penta(node,a1,ap1,ap2,am1,am2,q,ur1)      ! solution gets stored in ur1 - should always be u(r+1)
  if (iter == 1) then
    uo(1:node)=ur1(1:node)
  end if
  do i=1,node                         ! caclulate difference in solutions
    udiff=udiff + ABS(ur1(i)-ur(i))
  end do
  print *, udiff
  ur(1 : node) = ur1(1 : node)        ! store u(r+1) to u(r) for next step

! calculate the new load to add for next iteration step
! based on newest solution ur1 (u(r+1) 
! wherever it is <0, we need to add basin fill

  i=1 
  do while (i.le.node) 
    if (ur1(i) .lt. 0._QP) then
      q(i)=-1._QP*ur1(i)*g*pfill  ! fill value when the flexed plate is below zero - use ur1 as depth below zero and fill density
    else
      q(i)=0._QP                  ! sets remainder of q(i) to zero - because q(i) needs resetting because it was written over by the solution
    end if

    if (comp .eqv. .true.) then
      if (ur1(i).gt.0._QP) then
        q(i)=-1._QP*ur1(i)*g*pcrust  ! fill value when the flexed plate is below zero - use ur1 as depth below zero and fill density
!      else

!      q(i)=0._QP                  ! sets remainder of q(i) to zero - because q(i) needs resetting because it was written over by the solution
      end if
    end if 
    i=i+1
  end do

! add back in the original LOAD! 

  do i=1,node
    if (t(i).ne.0) then           ! checks places where there is load (t(i) ne 0)
      lf(i)=ur1(i)+abs(t(i))            ! lf is a variable to test if sum of original load and deflected plate is > or < 0 
      q(i)=t(i)*g*pload
      if (lf(i).lt.0) then         ! if lf < 0 then you need to add fill above the load whose top is < 0
        q(i)=q(i)-lf(i)*pfill*g      ! adds the additional fill above the load 
      end if 
     end if

!     if (tee(i).ne.h) then
!      q(i)=t(i)*g*pload
!     end if
  end do

! restore original boundary conditions

  q(1)=bc1
  q(2)=bc2
  q(node-1)=bc3
  q(node)=bc4

  if (bc3rd .eqv. .true.) then
    q(2)=0.5*bc2/vd(2)
    print *, 'q(2) = ', q(2) 
  end if         

  print *, iter, udiff

  iter=iter+1
end do

ur(1 : node) = ur1(1 : node)        ! store u(r+1) to u(r) for next step

199  continue

call plot(node,ur1,uo,vd,tee,t,q,lf,dx,pfill,pcrust,g,pint,code)

 50    continue     

 print *, MINVAL(ur1)
 print *, ur1(minloc(ur1)-1)

close(10)
close(11)
close(12)
close(13)
close(14)
close(15)
close(16)
close(17)
close(18)
close(19)
close(20)
close(21)

end

