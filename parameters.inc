!======================================================================
! For electroweak model:
! (following conventions in Achucarro and Vachaspati,
!  Phys. Rept. 327 (2000) 347-426.)
!
! nf=number of fields; nt=number of time steps;
! lat=semi-lattice size; dx=spatial spacing;
! alpha, beta etc. are parameters of the model
!
! Allows for different lattice sizes along different axes. 
!
! Fiducial lattice only goes up to lat-1. The last layer on the 
! lattice is a boundary layer where only some variables are calculated.
!======================================================================

      integer nf,latx,laty,latz,nt,nte,nsnaps,outerlayer,localsize
      integer nscalar,ngauge,nconstraints,nprocessx
      integer nf_snaps
      real*8 dx,gw,gy,lambda,vev,relaxparam, abs_tol
      real*8 pi
      real*8 Wangle
      integer senlen
!======================================================================
! For abelian-Higgs model: gy is the gauge coupling.
! Note: conventions like those in electroweak but truncated to hypercharge.
! model parameters:
!
! The number of fields (nf) includes real scalar fields, gauge fields,
! and the gauge constraints. E.g. for electroweak there are 4 real
! scalars, 4x4=16 gauge fields, and 4 gauge constraints (Lorenz gauge
! conditions. (gauge constraint equations are first order.)
! f(25) is assigned for |\phi|
!======================================================================

      parameter(pi=3.14159265359)
      parameter(nf=25,nscalar=4,ngauge=4*4,nconstraints=ngauge/4)

      parameter(gw=0.65,lambda=0.129,vev=1.0)
      parameter(Wangle=ASIN(sqrt(0.22)))
      parameter(gy=gw*tan(Wangle))

!======================================================================
! Simulation parameters:
!======================================================================
!Tolerance constraint on energy changes between subsequent iterations
      parameter(abs_tol=0.0001)

      parameter(latx=200,laty=latx,latz=latx,nprocessx=5,nt=100000,
     1 dx=0.25,localsize=(2*latx/nprocessx)+13,
     2 relaxparam=1.0)

!======================================================================
! Output parameter -- how many snapshots of the energy distribution?
!======================================================================

      parameter(nsnaps=nt/100)
      parameter(nf_snaps=nt/10)
!======================================================================
! Calculate energy in a sub-box with dimension -lat+outerlayer,+lat-outerlayer.
! outerlayer must be greater than equal to 1 (see energy subroutine).
!======================================================================

      parameter(outerlayer=2)

!      parameter(Wangle=0.488205)
!      parameter(Wangle=ATAN2(gy,gw))

!---Parameter that indicates the half length of slab in z-direction to
!---be integrated over for computing the string energy density
      parameter(senlen=20)
