!=======================================================================
! subroutine to compute and printout the magnetic fields
!=======================================================================

      subroutine electromagnetics(f,itime,pid,bb,fb,lb,rb,db,ub,
     1                   lxb,lyl,lzd,hatn)
            
      implicit none
      include 'parameters.inc'
      include 'initialParameters.inc'

      integer n,i,j,k,itime,pid,proci,procj,prock
      integer lxb,lyl,lzd,bb,fb,lb,rb,db,ub
      integer ilocal,jlocal,klocal,enbb,enfb,enlb,enrb,endb,enub
      integer lx, ly, lz
      integer dum_n
      real*8 x, y, z, dV
      real*8 f,fs
      real*8 dfdx,dfdy,dfdz,d2fdx2,d2fdy2,d2fdz2
      real*8 dfddx,dfddy,dfddz
      real*8 cd
      real*8 n_1,n_2,n_3
      real*8 g
      real*8 A,B
      real*8 cw,sw
      real*8 higgs_current
      complex*8 hatn

      dimension f(nf,1:localsize,1:localsize,1:localsize)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension dfddx(nf),dfddy(nf),dfddz(nf)
      dimension fs(4,0:3,0:3),d2fdx2(nf),d2fdy2(nf),d2fdz2(nf)
      dimension g(nf)
      dimension A(4,4)
      dimension B(3,1:localsize,1:localsize,1:localsize)
      dimension hatn(2,1:localsize,1:localsize,1:localsize)
!=======================================================================
! Covariant Derivatives:  cd(0:3,number of scalar fields):
!=======================================================================

      dimension cd(0:3,4)

      character (len=20) fname_B
      character (len=20) fname_B_xy

!       write(fname_B,"(A5,I3.3,A4)") "B_xz_",pid,".dat"
!       open(unit=1,file=trim(adjustl(fname_B)),position='append')
      
!       write(fname_B_xy,"(A5,I3.3,A4)") "B_xy_",pid,".dat"
!       open(unit=2,file=trim(adjustl(fname_B_xy)),position='append')

      cw = cos(Wangle)
      sw = sin(Wangle)

      !---Identifying the coordinates of the sublattice
!---Prefix l : sublattice coordinates
!---Prefix en : global lattice coordinates to constructed to omit
!---last lattice site 

      call processcoord(pid,proci,procj,prock)

       enbb=bb
       if(proci==nprocessx-1) then
              enfb=fb
       else
              enfb=fb-1
       endif

       enlb=lb
       if(procj==nprocessx-1) then
              enrb=rb
       else
              enrb=rb-1
       endif

       endb=db
       if(prock==nprocessx-1) then
              enub=ub
       else
              enub=ub-1
       endif
      

        do k=endb,enub
          do j=enlb,enrb
            do i=enbb,enfb
     
!     if one wants to only count the energy within a small region:
               ilocal = i+lxb-1
               jlocal = j+lyl-1
               klocal = k+lzd-1

               x=float(ilocal)*dx
               y=float(jlocal)*dx
               z=float(klocal)*dx

        do dum_n=1,nf
        g(dum_n)=f(dum_n,i,j,k)
        enddo


        call derivatives(f,i,j,k,dfdx,dfdy,dfdz,
     1                  d2fdx2,d2fdy2,d2fdz2,lxb,lyl,lzd)

!=======================================================================
! covariant derivative: cd(mu,a)=(D_\mu\Phi)^a (a=1 is the real 
! part of the upper component of the doublet; a=2 is the imag
! part; a=3 is real part of the lower component; a=4 is the
! imag part of the lower component.
!=======================================================================
         
        call covderiv(f,i,j,k,dfdx,dfdy,dfdz,cd)

!=======================================================================     
!find gauge field strengths (uses derivatives):
!=======================================================================

        call fieldstrength(f,i,j,k,dfdx,dfdy,dfdz,fs)

        n_1=-2.*(g(1)*g(3)+g(2)*g(4))/(g(25)**2)
        n_2=+2.*(g(2)*g(3)-g(1)*g(4))/(g(25)**2)
        n_3=(-g(1)**2-g(2)**2+g(3)**2+g(4)**2)/(g(25)**2)

!------- E&M gauge fields----------------------------------------------
!------- Using conventional definitions and indexes for EM tensor------!

!        higgs_current=

        A(4,3)=sw*(n_1*fs(1,3,2)+n_2*fs(2,3,2)+n_3*fs(3,3,2))
     1         +cw*fs(4,3,2)

        A(2,4)=sw*(n_1*fs(1,1,3)+n_2*fs(2,1,3)+n_3*fs(3,1,3))
     1         +cw*fs(4,1,3)

        A(3,2)=sw*(n_1*fs(1,2,1)+n_2*fs(2,2,1)+n_3*fs(3,2,1))
     1         +cw*fs(4,2,1)

        B(1,i,j,k)=A(4,3)
        B(2,i,j,k)=A(2,4)
        B(3,i,j,k)=A(3,2)


!-------printouts in xz plane
!         if(jlocal.eq.0) then
!         write(1,*) itime,ilocal,klocal,B(1,i,j,k),B(2,i,j,k),B(3,i,j,k)
!         endif
    
!         if(klocal.eq.0.or.klocal.eq.60.or.klocal.eq.(-60)) then
!         write(2,*) itime, ilocal, jlocal, klocal
!      1               ,B(1,i,j,k),B(2,i,j,k),B(3,i,j,k)
!         endif

                enddo
            enddo
        enddo

       !  close(1)
       !  close(2)

          return
      end
