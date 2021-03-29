MODULE plotflex

contains

subroutine plot(node,ur,uo,vd,tee,t,q,lf,dx,pfill,pcrust,g,pint,code)

implicit none

integer, parameter :: QP = selected_real_kind (32)

integer :: node,i,j,k,l,istore,check,pint
real (kind=QP) :: ur(node),uo(node),vd(node),tee(node),t(node),q(node),lf(node)
real (kind=QP) :: dx,pfill,pcrust,g,plotf
real (kind=QP), parameter :: eps=0.0001_QP
character*3 :: code

! OUTPUT FILES
if (code == 'ws3') then
  open (14, file='ws_udqp.txt')     ! output of solution (deflection of top surface) 
  open (15, file='ws_load.txt')     ! extra output for plotting
  open (16, file='ws_u-te.txt')     ! extra output for plotting  
  open (17, file='ws_fill.txt')     ! stores fill "hachuring"
  open (18, file='ws_uo.txt')       ! for writing the first solution before iteration starts
  open (19, file='ws_moment.txt')   ! moment along plate
  open (20, file='ws_logvd.txt')    !   
  open (21, file='ws_bending.txt')  ! bending along plate u''  
  open (22, file='ws_te-profile.txt') ! te for plotting   
  open (23, file='ws_fill1.txt')
end if

if (code == 'hs3') then
  open (14, file='hs_udqp.txt')     ! output of solution (deflection of top surface) 
  open (15, file='hs_load.txt')     ! extra output for plotting
  open (16, file='hs_u-te.txt')     ! extra output for plotting  
  open (17, file='hs_fill.txt')     ! stores fill "hachuring"
  open (18, file='hs_uo.txt')       ! for writing the first solution before iteration starts
  open (19, file='hs_moment.txt')   ! moment along plate
  open (20, file='hs_logvd.txt')    !   
  open (21, file='hs_bending.txt')  ! bending along plate u''  
  open (22, file='hs_te-profile.txt') ! te for plotting   
  open (23, file='hs_fill1.txt')
end if

plotf=dx/1000._QP
i=1
check=1
! plot files    ....  final iteration plot  
12 continue
  write (14,'(2(1pe15.5))') plotf*(i-1), ur(i)        ! top surface of lithosphere       u.txt
  if (i > 1 .and. i < node) then
    write (19,'(2(1pe15.5))') plotf*(i-1), abs(vd(i)*(ur(i-1)-2*ur(i)+ur(i+1)) / dx**2)
    write (21,'(2(1pe15.5))') plotf*(i-1), abs((ur(i-1)-2*ur(i)+ur(i+1)) / dx**2)
  else
    write (19,'(2(1pe15.5))') plotf*(i-1), 0._QP
    write (21,'(2(1pe15.5))') plotf*(i-1), 0._QP
  end if
  write (22,'(2(1pe15.5))') plotf*(i-1), tee(i)

! plot the load (file = 15)
       
  if (t(i).ne.0) then       ! use value of t(i) ne 0 to find load
    if (check == 1) then    ! first time load encountered (lhs of load) - t(i) ne 0 AND check = 1
      istore=i               ! store the node of lhs of current load
      write (15,'(a)') '>'  
      write (15,'(2(1pe15.5))') plotf*(i-1), ur(i)   ! so plot bottom of load on plate with value of ur(i) 
      check=2
    end if
    write (15,'(2(1pe15.5))') plotf*(i-1), ur(i)+t(i)   ! initial added load (not the fill)    load.txt
  end if
!  write (16,*) plotf*(i-1), ur(i)-h      ! bottom surface of lithosphere/crust ***as defined by ELASTIC thickness***  u-te.txt
  write (18,'(2(1pe15.5))') plotf*(i-1), uo(i)

  if (t(i).eq.0) then                       ! t(i) is applied load thickness. There should be infill or erosion only where there is no applied load. 
                                            ! This might actually be changed for additional modelling scenarios. Modifications would be made to flex-hs.f90 first, however.
   if (abs(q(i))>eps) then                  ! q(i) is basin fill/erosion generated force. If pcrust or pfill = 0, then q(i) = 0, hence this test prevents division by zero    
    if (ur(i).lt.0) then                    ! no load, basin fill, u -ve
      write (17,'(a)') '>'                                            ! fill hachuring 
      write (17,'(1pe15.5,i3)') plotf*(i-1), 0                        
      write (17,'(2(1pe15.5))') plotf*(i-1), -1._QP*q(i)/(pfill*g)    ! bottom of line given as load force (q) / (fill density*g)
      
    else                                    ! no load, erosion, u +ve
      
      write (23,'(a)') '>'                                            ! fill hachuring
      write (23,'(1pe15.5,i3)') plotf*(i-1), 0                        ! bottom of line
      write (23,'(2(1pe15.5))') plotf*(i-1), -1._QP*q(i)/(pcrust*g)   ! top of line given as load force (q) / (crust density*g)
    end if
   end if 
    if (check == 2) then    ! complete the edge of the load - where t(i) eq 0 for first time (check = 2) means you've reached the rhs of the load
      write (15,'(2(1pe15.5))') plotf*(i-2), ur(i-1) ! so plot bottom corner of load at top of plate at node i-1 !!!
      do j=i-2,istore+1,-1    ! write in ur from rhs back to lhs of current load (nodes i-2 to istore+1)
        write (15,'(2(1pe15.5))') plotf*(j-1), ur(j)
      end do
      check=1
    end if

  end if
  if (lf(i).lt.0) then
    write (17,'(a)') '>'
    write (17,'(1pe15.5,i3)') plotf*(i-1), 0               ! fill hachuring
    write (17,'(2(1pe15.5))') plotf*(i-1), lf(i)    ! fill hachuring  where there is a load that is below zero
  end if 
  i=i+pint
  if (i > node) then
    goto 50
  else
    goto 12
  end if

 50 continue    
 
close(14)
close(15)
close(16)
close(17)
close(18)
close(19)
close(20)
close(21)
close(22)
close(23)

end subroutine plot

end module plotflex       
