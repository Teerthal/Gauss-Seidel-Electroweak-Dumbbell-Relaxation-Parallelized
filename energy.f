!=======================================================================
! SU(2)xU(1) energy:
!=======================================================================

!=======================================================================
!This subroutine computes energies in the all field components of the theory
!and also performs volume integrations within the parallelized 
!domains. The total energy of each field component within each domain is  
!passed in the energy array along with the iteration stamp
!Includes commented out example routines to save energy densities in files
!=======================================================================

      subroutine uenergy(f,energy,itime,pid,lxb,lyl,lzd,
     1                          bb,fb,lb,rb,db,ub)
      implicit none
      include 'parameters.inc'
      include 'initialParameters.inc'

      integer n,i,j,k,itime,pid,proci,procj,prock
      integer lxb,lyl,lzd,bb,fb,lb,rb,db,ub
      integer ilocal,jlocal,klocal,enbb,enfb,enlb,enrb,endb,enub
      integer lx, ly, lz
      integer dum_n
      real*8 x, y, z, dV
      real*8 kineticEnergyPhi, gradientEnergyPhi, electricEnergyW
      real*8 magneticEnergyW, electricEnergyY, magneticEnergyY
      real*8 potentialEnergy, energyDensity
      real*8 totalKEPhi, totalGEPhi, totalEEW, totalMEW, totalEEY
      real*8 totalMEY, totalPE, totalEnergy
      real*8 totalstringenergy
      real*8 f
      real*8 dfdx,dfdy,dfdz,fs,d2fdx2,d2fdy2,d2fdz2
      real*8 dfddx,dfddy,dfddz
      real*8 cd
      real*8 energyInSphere,r, energy(13), energyBoundary

      dimension f(nf,1:localsize,1:localsize,1:localsize)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension dfddx(nf),dfddy(nf),dfddz(nf)
      dimension fs(4,0:3,0:3),d2fdx2(nf),d2fdy2(nf),d2fdz2(nf)

      character (len=20) fname_phi_xy
      character (len=20) fname_phi_xz
!=======================================================================
! Covariant Derivatives:  cd(0:3,number of scalar fields):
!=======================================================================

      dimension cd(0:3,4)

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
      
        
!=======================================================================
! total scalar KE, gradient energy; total gauge electric, magnetic; 
! total scalar potential; total energy.
!=======================================================================

        totalKEPhi=0.
        totalGEPhi=0.
        totalEEW=0.
        totalMEW=0.
        totalEEY=0.
        totalMEY=0.
        totalPE=0.
        totalEnergy=0.
        energyInSphere=0.
        energyBoundary=0.
        totalstringenergy=0.


       write(fname_phi_xy,"(A7,I3.3,A4)") "phi_xy_",pid,".dat"
       open(unit=53,file=trim(adjustl(fname_phi_xy)),position='append')
  
       write(fname_phi_xz,"(A7,I3.3,A4)") "phi_xz_",pid,".dat"
       open(unit=54,file=trim(adjustl(fname_phi_xz)),position='append')


!=======================================================================
!     volume element:
!=======================================================================

        dV=dx**3
       
        do k=endb,enub
          do j=enlb,enrb
            do i=enbb,enfb
     
               ilocal = i+lxb-1
               jlocal = j+lyl-1
               klocal = k+lzd-1

               x=float(ilocal)*dx
               y=float(jlocal)*dx
               z=float(klocal)*dx

        call derivatives(f,i,j,k,dfdx,dfdy,dfdz,
     1                  d2fdx2,d2fdy2,d2fdz2,lxb,lyl,lzd,itime)

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

        kineticEnergyPhi= (cd(0,1)**2+cd(0,2)**2
     1                         +cd(0,3)**2+cd(0,4)**2)

        gradientEnergyPhi= (cd(1,1)**2+cd(2,1)**2+cd(3,1)**2
     1                       +cd(1,2)**2+cd(2,2)**2+cd(3,2)**2
     2                       +cd(1,3)**2+cd(2,3)**2+cd(3,3)**2
     3                       +cd(1,4)**2+cd(2,4)**2+cd(3,4)**2)

        electricEnergyW=0.5*(fs(1,0,1)**2+fs(1,0,2)**2+fs(1,0,3)**2
     1                    +fs(2,0,1)**2+fs(2,0,2)**2+fs(2,0,3)**2
     1                    +fs(3,0,1)**2+fs(3,0,2)**2+fs(3,0,3)**2)

        magneticEnergyW=0.5*(fs(1,1,2)**2+fs(1,1,3)**2+fs(1,2,3)**2
     1                    +fs(2,1,2)**2+fs(2,1,3)**2+fs(2,2,3)**2
     2                    +fs(3,1,2)**2+fs(3,1,3)**2+fs(3,2,3)**2)

!Check Y energies
        electricEnergyY=0.5*(fs(4,0,1)**2+fs(4,0,2)**2+fs(4,0,3)**2)
        
        magneticEnergyY=0.5*(fs(4,1,2)**2+fs(4,1,3)**2+fs(4,2,3)**2)

        potentialEnergy=lambda*(f(1,i,j,k)**2+f(2,i,j,k)**2
     1                       +f(3,i,j,k)**2+f(4,i,j,k)**2-vev**2)**2

! ======================================================================
     
       energyDensity=kineticEnergyPhi+gradientEnergyPhi+electricEnergyW+
     1   magneticEnergyW+electricEnergyY+magneticEnergyY+potentialEnergy
        
        totalKEPhi=totalKEPhi+kineticEnergyPhi*dV
        totalGEPhi=totalGEPhi+gradientEnergyPhi*dV
        totalEEW=totalEEW+electricEnergyW*dV
        totalMEW=totalMEW+magneticEnergyW*dV
        totalEEY=totalEEY+electricEnergyY*dV
        totalMEY=totalMEY+magneticEnergyY*dV
        totalPE=totalPE+potentialEnergy*dV
        totalEnergy=totalEnergy+energyDensity*dV

!=======================================================================
! energy within sphere -- useful for single monopole:
!        r=sqrt((x-xm)**2+(y-ym)**2+(z-zm)**2)
! energy within sphere -- for mmbar the center is at origin:
!=======================================================================

        r=sqrt(x**2+y**2+z**2)

        if(r.lt.float(latx-1)*dx) then
        energyInSphere=energyInSphere+energyDensity*dV
        endif

        if((abs(ilocal).ge.(latx)).or.(abs(jlocal).ge.(latx))
     1              .and.(abs(klocal).ge.(latx))) then
        energyBoundary=energyBoundary+energyDensity*dV
        endif

!=======================================================================
!Energy within slab of length senlen
!=======================================================================
        
        if(abs(klocal).lt.(senlen)) then
        totalstringenergy=totalstringenergy+energyDensity*dV
        endif        

        if(jlocal.eq.0) then
         write(54,"(I5,I5,I5,E25.16)")
     -   itime,ilocal,klocal,f(25,i,j,k)
        endif

        if(klocal.eq.0.or.klocal.eq.60.or.klocal.eq.(-60)) then
          write(53,"(I5,I5,I5,E25.16)")
     -    itime,ilocal,jlocal,f(25,i,j,k)
        endif

                enddo
            enddo
        enddo

        close(53)
        close(54)
        
        energy(1)=itime
        energy(2)=totalKEPhi
        energy(3)=totalGEPhi
        energy(4)=totalEEW
        energy(5)=totalMEW
        energy(6)=totalEEY
        energy(7)=totalMEY
        energy(8)=totalPE
        energy(9)=totalEnergy
        energy(10)=float(latx-1)*dx
        energy(11)=energyInSphere
        energy(12)=energyBoundary
        energy(13)=totalstringenergy

          return
      end
