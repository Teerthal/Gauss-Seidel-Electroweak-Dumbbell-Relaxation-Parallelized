!=======================================================================
! Relax monopole-antimonopole field configuration using Gauss-Seidel 
! method.
!=======================================================================

! begin{declarations}

      use mpi
      implicit none

!=======================================================================
! begin{input}
! model-dependent parameters:

      include 'parameters.inc'

! problem-dependent parameters:

      include 'initialParameters.inc'

! end{input}
! ======================================================================

!     Declare all fields and time derivatives -- nf is the number of
!     fields; f denotes field. 

      real*8 f 
!      real*8 gauge, mf
      real*8 energyinitial, totalEnergy
      real*8 energy(13), energy0(13)
      real*8 energymin
      complex*8 hatn
      complex*8 dxhatn, dyhatn, dzhatn
      complex*8 d2xhatn, d2yhatn, d2zhatn
      real*8 csNumber

      dimension f(nf,1:localsize,1:localsize,1:localsize)
!      dimension gauge(4,1:localsize,1:localsize,1:localsize)
!      dimension mf(3,1:localsize,1:localsize,1:localsize)
      dimension hatn(2,1:localsize,1:localsize,1:localsize)
      dimension dxhatn(2,1:localsize,1:localsize,1:localsize)
      dimension dyhatn(2,1:localsize,1:localsize,1:localsize)
      dimension dzhatn(2,1:localsize,1:localsize,1:localsize)
      dimension d2xhatn(2,1:localsize,1:localsize,1:localsize)
      dimension d2yhatn(2,1:localsize,1:localsize,1:localsize)
      dimension d2zhatn(2,1:localsize,1:localsize,1:localsize)

      integer it, itime,i,j,k,l,ierr,np,pid,interval
      integer lxb,lyl,lzd,bb,fb,lb,rb,db,ub
      integer proci,procj,prock
      integer n
      integer ilocal,jlocal,klocal

      character (len=40) fname
      character (len=40) parameter_file
      character (len=40) d_name


! end{declarations}

        call mpi_init(ierr)
        call mpi_comm_rank(mpi_comm_world, pid, ierr)
        call mpi_comm_size(mpi_comm_world, np, ierr)

        call boundindices(pid,lxb,lyl,lzd,bb,fb,lb,rb,db,ub)
!        print*,pid,lxb,lyl,lzd,bb,fb,lb,rb

        if(pid==0) then
        print *, ' no. of fields =                       ', nf
        print *, ' lattice size =                        ', latx
        print *, ' lattice spacing =                     ', dx
        print *, ' Gauss-Seidel relaxparam =             ', relaxparam
        print *, ' number of snap shots =                ', nsnaps
        print *, ' gauge function parameter =            ', gp2
        print *, ' theory parameters: gw,gy,lambda,vev = ', gw, gy,
     1                  lambda,vev
        print *, 'number of processors = ', np
        print *, ' '
        endif

        if(pid==0) then
        write(parameter_file,"(A16)") "export_paras.txt"
        open(unit=99,file=trim(adjustl(parameter_file))
     1      ,position='append',status='replace')
        write(99,*) latx
        write(99,*) nprocessx
        write(99,*) nte
        write(99,*) nsnaps_e
        close(99)
        endif

!=======================================================================
! begin{initial conditions}
!=======================================================================

        call initialconditions(f,pid,
     1                          lxb,lyl,lzd,bb,fb,lb,rb,db,ub,np,
     2                          hatn,dxhatn,dyhatn,dzhatn,
     3                          d2xhatn,d2yhatn,d2zhatn)
     
!=======================================================================
! end{initial conditions}
!=======================================================================
        
!=======================================================================
! begin{write out initial energies and chern simons number}
! 'it' in subroutine energy is the time step (zero right now):
!=======================================================================

        it=0
        call uenergy(f,energy,it,pid,lxb,lyl,lzd,bb,fb,lb,rb,db,ub)

!        call chernsimons(f,csNumber,it,pid,lxb,lyl,lzd,
!     1                          bb,fb,lb,rb,db,ub)

!=======================================================================
! remember initial energy and chern simons number in box:
!=======================================================================
        call printenergy(energy,pid,np)
!        call printCS(csNumber,pid,np,it)
        energyinitial=energy(9)

!=======================================================================
! initalize energymin (used to stop after the energy has become minimum)
!=======================================================================

         energymin=energyinitial

!=======================================================================
! Printout initial magnetic field
!=======================================================================
!        call electromagnetics(f,it,pid,bb,fb,lb,rb,db,ub,
!     1                   lxb,lyl,lzd,hatn)

!=======================================================================
! begin{evolution}
!=======================================================================

        do itime=1,nt
        
!=======================================================================
!Subrountines to call relaxation subroutine as well as as mpi send and
!receive relevant ghost cell data
!=======================================================================

        call evolveeuler(f,pid,np,hatn,dxhatn,dyhatn,dzhatn,
     1                 d2xhatn,d2yhatn,d2zhatn,
     2                 lxb,lyl,lzd,bb,fb,lb,rb,db,ub,itime)
        

!Compute and print energies at specified intervals
        interval=int(nt/nsnaps)

        if(int(nt/nsnaps).eq.0) interval=1
        if(mod(itime,interval).eq.0.or.itime.eq.nt) then

        call uenergy(f,energy,itime,pid,lxb,lyl,lzd,bb,fb,lb,rb,db,ub)

        call printenergy(energy,pid,np)

!        call electromagnetics(f,itime,pid,bb,fb,lb,rb,db,ub,
!     1                   lxb,lyl,lzd,hatn)


        totalEnergy=energy(9)

!       print*, '---------------------------'
!       print*, 'total energy', totalEnergy
!       print*, '---------------------------'

!=======================================================================
! Stop if energy minimum has been reached. Once the energy
! starts growing, it is seen to grow rapidly, signaling an
! instability. So best to stop if energy grows:
!=======================================================================
        if(pid==0) then
        print*,'------------------------'
        print*,'Subsequent energy difference'
        print*,abs(totalEnergy-energymin)/energymin
        print*,'------------------------'
        endif
        
        if(totalEnergy.le.energymin.and.
     1          abs(totalEnergy-energymin)/energymin.gt.abs_tol) then

                energymin=totalEnergy
!Dump data arrays evey energy compute step
!and replace previous until reaching the 2nd to 
!last iteration

        write(fname,'(a,I3.3,A4)') "rdata/fin_f_",pid,".dat"
        open(unit=pid,file=trim(adjustl(fname)),form='unformatted')
        write(pid) f
        close(pid)

!Could uncomment and run the following statements to check if the file
!saves were successful
!        open(unit=pid,file=fname,form='unformatted',
!      1       action="read")

!         read(pid) f
!         close(pid)
        
        else

        if(pid==0) then
                print*, ' totalEnergy.gt.energymin =  STOP '
        endif
                exit
        endif

        endif

        enddo

      end
