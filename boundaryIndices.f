!------Subrountines for the boundary indices------
! Here the processcoord subroutine converts the pid to
! cartesian processor coordinates i.e it assigns a sublattice
! to each processor. The subroutine itself assigns global lattice coords
! corresponding to the local coordinates using processor coordinates (procii,..)
!Local boundary coordinates are shifted by 7 for ghost cell space


        subroutine boundindices(pid,lxb,lyl,lzd,bb,fb,lb,rb,db,ub)
                implicit none

        include 'parameters.inc'
        integer lxb,lxf,lyl,lyr,lzd,lzu
        integer bb,fb,lb,rb,db,ub
        integer pid,proci,procj,prock

        call processcoord(pid,proci,procj,prock)

!Dictionary
! Prefixes l : global integer cartesian coordinates (x,y,z)
! b:boundary
! Prefixes before b(boundary) :
!--------------------------------------------|
!Cartesian coordinate  |   boundary sides    |
!--------------------------------------------|
!          x           |  b:back  ; f:front ;|
!          y           |  l:left  ; r:right ;|
!          z           |  d:down  ; u:up     |
!--------------------------------------------

!-------------------x-indices----------------------------

                lxb=-latx+proci*((2*latx)/nprocessx)-6
                lxf=-latx+(proci+1)*((2*latx)/nprocessx)+6
                bb=7
                fb=((2*latx)/nprocessx)+7

!-------------------y-indices----------------------------

                lyl=-latx+procj*((2*latx)/nprocessx)-6
                lyr=-latx+(procj+1)*((2*latx)/nprocessx)+6
                lb=7
                rb=((2*latx)/nprocessx)+7

!-------------------z-indices----------------------------

                lzd=-latx+prock*((2*latx)/nprocessx)-6
                lzu=-latx+(prock+1)*((2*latx)/nprocessx)+6
                db=7
                ub=((2*latx)/nprocessx)+7

!---Testing exectuion----Checked for consistency----
!        print*,' ======================= '
!                print*,pid,bb,lxb,fb,lxf
!        print*,' ======================= '
!------------------------------------------

!---The ones below seem phased out-----
!-------------------x-indices----------------------------

!       if(proci==0) then
!
!               lxb=-latx
!               lxf=-latx+(2*latx)/nprocessx+13-1
!               bb=1
!               fb=((2*latx)/nprocessx)+1
!
!       elseif(proci==nprocessx-1) then
!
!               lxb=latx-((2*latx)/nprocessx)-13+1
!               lxf=latx
!               bb=13
!               fb=((2*latx)/nprocessx)+13
!
!       else
!
!               lxb=-latx+proci*((2*latx)/nprocessx)-6
!               lxf=-latx+(proci+1)*((2*latx)/nprocessx)+6
!               bb=7
!               fb=((2*latx)/nprocessx)+7
!
!       endif
!
!-------------------y-indices----------------------------
!
!       if(procj==0) then
!
!               lyl=-latx
!               lyr=-latx+(2*latx)/nprocessx+13-1
!               lb=1
!               rb=((2*latx)/nprocessx)+1
!
!       elseif(procj==nprocessx-1) then
!
!               lyl=latx-((2*latx)/nprocessx)-13+1
!               lyr=latx
!               lb=13
!               rb=((2*latx)/nprocessx)+13
!
!       else

!               lyl=-latx+procj*((2*latx)/nprocessx)-6
!               lyr=-latx+(procj+1)*((2*latx)/nprocessx)+6
!               lb=7
!               rb=((2*latx)/nprocessx)+7
!
!       endif

!-------------------z-indices----------------------------
!
!       if(prock==0) then
!
!               lzd=-latx
!               lzu=-latx+(2*latx)/nprocessx+13-1
!               db=1
!               ub=((2*latx)/nprocessx)+1
!
!       elseif(prock==nprocessx-1) then
!
!               lzd=latx-((2*latx)/nprocessx)-13+1
!               lzu=latx
!               db=13
!               ub=((2*latx)/nprocessx)+13
!
!       else

!               lzd=-latx+prock*((2*latx)/nprocessx)-6
!               lzu=-latx+(prock+1)*((2*latx)/nprocessx)+6
!               db=7
!               ub=((2*latx)/nprocessx)+7
!
!       endif


        return
        end
