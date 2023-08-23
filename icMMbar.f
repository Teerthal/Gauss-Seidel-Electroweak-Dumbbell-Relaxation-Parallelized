!=======================================================================
! Initial conditions for the dumbbell in su2u1 model. 
! See arXiv:2302.04886 for explanation for the 
! profile makeup and choices
!=======================================================================

!----------Variable Definitions---------------------------
!--f = field array 
!--localsize = local array size
!--indices run from 1:localsize
!--pprofilem: Radial Monopole profile function
!--pprofilea: Radial Anti-Monopole profile function
!--sprofile: String profile function
      subroutine initialconditions(f,pid,
     1                            lxb,lyl,lzd,bb,fb,lb,rb,db,ub,np,
     2                            hatn,dxhatn,dyhatn,dzhatn,
     3                          d2xhatn,d2yhatn,d2zhatn)

      implicit none
      include 'parameters.inc'
      include 'initialParameters.inc'
      integer i,j,k,n
      integer izm
      integer dum_n
      integer spi
      real*8 itw
      real*8 xa,ya,za
      real*8 x,y,z,xgm,ygm,zgm,xga,yga,zga,rm,ra,rhom,rhoa
      real*8 sm,sa,ms,mv
      real*8 pprofilem,pprofilea,wprofilem,wprofilea,sprofile
      real*8 yprofilem, yprofilea
      real*8 f
!      real*8 cf
      real*8 phisq
      real*8 twist,ctw,stw,ctw2,stw2
      real*8 correctionm,correctiona
      real*8 theta_m,theta_a, sph_phi
      real*8 jac, jac_2
      real*8 xua,yua,zua,xum,yum,zum,rum,rua
      real*8 cy_rm, cy_ra
      real*8 cw,sw, S2,S2b,C2,C2b
      integer bx,fx,ly,ry,dz,uz,lxb,lyl,lzd
      integer iglobal,jglobal,kglobal,np,pid,bb,fb,lb,rb,db,ub
      complex*8 hatn
      complex*8 dxhatn, dyhatn, dzhatn
      complex*8 d2xhatn, d2yhatn, d2zhatn
      complex*8 db_phi_x,db_phi_y,db_phi_z
!      real*8 dfdx,dfdy,dfdz,fs,d2fdx2,d2fdy2,d2fdz2
!      real*8 dfddx,dfddy,dfddz
      real*8 n_1, n_2, n_3
      real*8 n_1_x,n_1_y,n_1_z
      real*8 n_2_x,n_2_y,n_2_z
      real*8 n_3_x,n_3_y,n_3_z
      real*8 nx_mag,ny_mag,nz_mag
      real*8 tem_z
      dimension jac(2,3,3)
      dimension jac_2(2,3,3)
      dimension f(nf,1:localsize,1:localsize,1:localsize)
!      dimension cf(nf,1:localsize,1:localsize,1:localsize)
!      dimension dfdx(nf),dfdy(nf),dfdz(nf)
!      dimension dfddx(nf),dfddy(nf),dfddz(nf)
      dimension hatn(2,1:localsize,1:localsize,1:localsize)
      dimension dxhatn(2,1:localsize,1:localsize,1:localsize)
      dimension dyhatn(2,1:localsize,1:localsize,1:localsize)
      dimension dzhatn(2,1:localsize,1:localsize,1:localsize)
      dimension d2xhatn(2,1:localsize,1:localsize,1:localsize)
      dimension d2yhatn(2,1:localsize,1:localsize,1:localsize)
      dimension d2zhatn(2,1:localsize,1:localsize,1:localsize)
!=======================================================================
! This section is for running the code on several nodes at the
! same time, each code with a different value of the parameters.
! read scan parameters and create files to write out successful runs 
!=======================================================================

        character*500 filename2
        character (len=200) winfile
        character (len=200) itwist
        character (len=200) izmono
        character (len=20) filename3
        character (len=20) fnamei
        character (len=20) fnamehatn

!        itw = 0.0 
!        izm = 60
        call get_command_argument(1,winfile)
        read(winfile,'(a500)')filename2

        call get_command_argument(2,itwist)
        read(itwist,*) itw
        twist=(itw)*3.1416/6.

        call get_command_argument(3,izmono)
        read(izmono,*) izm
        zm=(float(izm)+0.5)*dx
     


        if(pid==0) then
        print*,' ================= '
        print*,' lambda, twist, zm ', lambda,twist,zm
        print*,' ================= '
        endif

!=======================================================================
!       initialize:
!=======================================================================

        do k=1,localsize
            do j=1,localsize
                do i=1,localsize

                do n=1,nf
                  f(n,i,j,k)=0.
                enddo
                enddo
            enddo
        enddo



!=======================================================================
! fields for SU(1)xU(1):
!     f(1)=phi^1,f(2)=phi^2,f(3)=phi^3,f(4)=phi^4
!     f(5)=W^1_0, f(6)=W^1_1, f(7)=W^1_2, f(8)=W^1_3,
!     f(9)=W^2_0, f(10)=W^2_1, f(11)=W^2_2, f(12)=W^2_3,
!     f(13)=W^3_0, f(14)=W^3_1, f(15)=W^3_2, f(16)=W^3_3,
!     f(17)=B_0, f(18)=B_1, f(19)=B_2, f(20)=B_3 
!     gauge fixing functions (Num Rel): f(21)=\partial_i W^1_i, 
!     f(22)=\partial_i W^2_i, f(23)=\partial_i W^3_i, 
!     f(24)=\partial_iB_i
!     time derivatives:
!     fd(i)=time derivative of f(i) for i=1-4 (scalar fields),
!     fd(i)=electric fields(=W^a_{0j}) for i=5-20 (gauge fields).
!     fd(24)=time derivative of |phi|
!=======================================================================

!=======================================================================
! scalar mass in model:
!=======================================================================

        ms=sqrt(2.*lambda)*vev
        mv=gw*vev

        if(pid==0) then
        print*, ' ms, mv, & ms/mv: ', ms,mv,ms/mv
        print*,' ======================= '
        print*, ' '
        endif

!=======================================================================
! Initial location of monopole is specified in initialParameters.inc
! (or passed on as an argument when executing) and is chosen so that 
! it is not on a lattice site.
! Antimonopole will be at (xm,ym,-zm):
!=======================================================================

        xa=xm
        ya=ym
        za=-zm	

! cosine and sine of twist:

        ctw=cos(twist)
        stw=sin(twist)
        ctw2=cos(twist/2.)
        stw2=sin(twist/2.)

        cw = cos(Wangle)
        sw = sin(Wangle)

!=======================================================================
! start: loop over larger lattice (+1) for scalar fields:
!=======================================================================

!        write(filename3,"(A4,I3.3,A4)") "ini_",pid,".dat"
!        open(unit=10,file=trim(adjustl(filename3)))

        do k=1,localsize
            do j=1,localsize
                do i=1,localsize

! ------------------ begin: setup variables -----------------
!               (x,y,z) for lattice site:
!--local labels are actually global coordinates
!--x,y,z are global coordinates
!--xgm, ... xga : global coordinates for monopole(m) and antimonopole(a) 
!---Boosts are in the x-direction
!--h are the higgs isovector direction as in eq(12) of 1705.0309

                iglobal = i+lxb-1
                jglobal = j+lyl-1
                kglobal = k+lzd-1

                x=float(iglobal)*dx
                y=float(jglobal)*dx
                z=float(kglobal)*dx

! Unboosted coordinates

                xum = (x-xm)
                yum = (y-ym)
                zum = (z-zm)
                xua = (x-xa)
                yua = (y-ya)
                zua = (z-za)
                rum = sqrt(xum**2+yum**2+zum**2)
                rua = sqrt(xua**2+yua**2+zua**2)
                cy_rm = sqrt(xum**2+yum**2)
                cy_ra = sqrt(xua**2+yua**2)

! Spherical coordinates in terms of non-boosted
! cartesian coordinates                

                theta_m = acos(zum/rum)
                theta_a = acos(zua/rua)
! Since y and x coordinates are the same for the monpole and
! antimonopole, we use monopole-centered coordinates here to
! define the azimuthal angle. 

                sph_phi = atan2(yum,xum)

! dimensionless distances (in units of vector mass):

                sm=vev*rum
                sa=vev*rua

! -- end: setup variables ------
!
!--begin scalar profile -----
! monopole and antimonopole scalar profile functions:
                if(sm.ne.0.) then
                pprofilem=vev*(1./tanh(2.*sm)-
     1           (1.+2.*ms*sm)*exp(-2.*ms*sm)/(2.*sm))
!                pprofilem=(1.-exp(-sm))
                else
                pprofilem=0.
                endif

                if(sa.ne.0.) then
                pprofilea=vev*(1./tanh(2.*sa)-
     1           (1.+2.*ms*sa)*exp(-2.*ms*sa)/(2.*sa))
!                pprofilem=(1.-exp(-sa))
                else
                pprofilea=0.
                endif

                sprofile=(1.-1.*((0.5*(tanh(-z+zm+szoff))
     1                  +0.5*(tanh(z-za+szoff)))/(tanh(zm+szoff)))
     2                  *(((0.5*tanh(-cy_ra+dx)
     1                   +0.5*tanh(cy_ra+dx))/tanh(dx))))

!                  sprofile=(1.-(cy_ra/sinh(cy_ra)))

!                sprofile=((-(((tanh(-z+zm+0.5*dx)+1.)/2.
!     1                    +(tanh(z-za+0.5*dx)+1.)/2.-1.)/(
!     2                      tanh(zm+0.5*dx)))
!     3                   *((tanh(-cy_rm+str_rad)+1.)/2.+
!     4                    (tanh(cy_rm+str_rad)+1.)/2.-1.)/(
!     5                     tanh(str_rad)))+1.)


!                sprofile=(1.-((0.5*tanh(-cy_ra+dx)
!     1                   +0.5*tanh(cy_ra+dx))/tanh(dx)))

!---------------------end scalar profiles and derivatives --------------

! Profiles for the scalar
! The if statements make sure that no Higgs zeros are 
! on the lattice size. This is unncessary for this setup
! since the initial profile and Higgs ansatz ensures that no zeros 
! are evaluated at a lattice size. This requires the
! initial MMbar location parameters in the initialParamaters.inc file
! to be appropriately chosen

        if((rum*rua).ne.0.) then
    
!          f(1,i,j,k)=sprofile*pprofilem*pprofilea*(
!     1                 sin(theta_m/2.)*sin(theta_a/2.)*ctw
!    2                +cos(theta_m/2.)*cos(theta_a/2.))/(vev)

!          f(2,i,j,k)=sprofile*pprofilem*pprofilea*(
!     1                  sin(theta_m/2.)*sin(theta_a/2.)*
!     2                  stw)/vev

!          f(3,i,j,k)=sprofile*pprofilem*pprofilea*(
!     1                  sin(theta_m/2.)*cos(theta_a/2.)*cos(sph_phi)
!     2                  -cos(theta_m/2.)*sin(theta_a/2.)
!     3                  *cos(sph_phi-twist))/vev

!          f(4,i,j,k)=sprofile*pprofilem*pprofilea*(
!     1                  sin(theta_m/2.)*cos(theta_a/2.)*sin(sph_phi)
!     2                  -cos(theta_m/2.)*sin(theta_a/2.)
!     3                  *sin(sph_phi-twist))/vev

          hatn(1,i,j,k)=cmplx(
     1                 sin(theta_m/2.)*sin(theta_a/2.)*ctw
     1                +cos(theta_m/2.)*cos(theta_a/2.),
     2                  (sin(theta_m/2.)*sin(theta_a/2.)*stw))

          hatn(2,i,j,k)=cmplx(
     1       sin(theta_m/2.)*cos(theta_a/2.)*cos(sph_phi)
     1      -cos(theta_m/2.)*sin(theta_a/2.)*cos(sph_phi-twist),
     3       sin(theta_m/2.)*cos(theta_a/2.)*sin(sph_phi)
     4      -cos(theta_m/2.)*sin(theta_a/2.)*sin(sph_phi-twist))

        f(1,i,j,k)=sprofile*pprofilea*pprofilem*REAL(hatn(1,i,j,k))
     1            /vev
        f(2,i,j,k)=sprofile*pprofilea*pprofilem*AIMAG(hatn(1,i,j,k))
     1            /vev
        f(3,i,j,k)=sprofile*pprofilea*pprofilem*REAL(hatn(2,i,j,k))
     1            /vev
        f(4,i,j,k)=sprofile*pprofilea*pprofilem*AIMAG(hatn(2,i,j,k))
     1            /vev
        f(25,i,j,k)=sprofile*pprofilem*pprofilea/vev

        else
        
        f(1,i,j,k)=0.
        f(2,i,j,k)=0.
        f(3,i,j,k)=0.
        f(4,i,j,k)=0.
        f(25,i,j,k)=0.

!-------Place holders----------!
        hatn(1,i,j,k)=cmplx(0.,0.)
        hatn(2,i,j,k)=cmplx(0.,0.)

        endif
                 enddo
             enddo
         enddo


!--------derivatives of hatn------------!


        do k=1,localsize
            do j=1,localsize
                do i=1,localsize

! ------------------ begin: setup variables -----------------
!               (x,y,z) for lattice site:
!--local labels are actually global coordinates
!--x,y,z are global coordinates
!--xgm, ... xga : global coordinates for monopole(m) and antimonopole(a) 
!---Boosts are in the x-direction
!--h are the higgs isovector direction as in eq(12) of 1705.0309

                iglobal = i+lxb-1
                jglobal = j+lyl-1
                kglobal = k+lzd-1

                x=float(iglobal)*dx
                y=float(jglobal)*dx
                z=float(kglobal)*dx

! Unboosted coordinates

                xum = (x-xm)
                yum = (y-ym)
                zum = (z-zm)
                xua = (x-xa)
                yua = (y-ya)
                zua = (z-za)
                rum = sqrt(xum**2+yum**2+zum**2)
                rua = sqrt(xua**2+yua**2+zua**2)
                cy_rm = sqrt(xum**2+yum**2)
                cy_ra = sqrt(xua**2+yua**2)

! Spherical coordinates in terms of non-boosted
! cartesian coordinates                

                theta_m = acos(zum/rum)
                theta_a = acos(zua/rua)
! Since y and x coordinates are the same for the monpole and
! antimonopole, we use monopole-centered coordinates here to
! define the azimuthal angle.

                sph_phi = atan2(yum,xum)

! dimensionless distances (in units of vector mass):

                sm=vev*rum
                sa=vev*rua

! Jacobian elements
! First index:1/2 -> monopole/antimonopole centered coordinates
! These are the elements corresponding to the partial derivatives of
! the spherical coordinates w.r.t the un-boosted cartesian coordinates.


                jac(1,1,1) = xum/rum
                jac(1,1,2) = yum/rum
                jac(1,1,3) = zum/rum
                jac(1,2,1) = xum*zum/(rum**2*sqrt(xum**2+yum**2))
                jac(1,2,2) = yum*zum/(rum**2*sqrt(xum**2+yum**2))
                jac(1,2,3) = -sqrt(xum**2+yum**2)/(rum**2)
                jac(1,3,1) = -yum/(xum**2+yum**2)
                jac(1,3,2) = xum/(xum**2+yum**2)
                jac(1,3,3) = 0.

                jac(2,1,1) = xua/rua
                jac(2,1,2) = yua/rua
                jac(2,1,3) = zua/rua
                jac(2,2,1) = xua*zua/(rua**2*sqrt(xua**2+yua**2))
                jac(2,2,2) = yua*zua/(rua**2*sqrt(xua**2+yua**2))
                jac(2,2,3) = -sqrt(xua**2+yua**2)/(rua**2)
                jac(2,3,1) = -yua/(xua**2+yua**2)
                jac(2,3,2) = xua/(xua**2+yua**2)
                jac(2,3,3) = 0.

                
                jac_2(1,1,1) = (yum**2+zum**2)/(rum**3)
                jac_2(1,1,2) = (xum**2+zum**2)/(rum**3)
                jac_2(1,1,3) = (xum**2+yum**2)/(rum**3)
                jac_2(1,2,1) = (zum*((-2.*xum**4)-((xum**2)*(yum**2))
     1                         +(yum**4)+((yum**2)*(zum**2))))
     2                         /((cy_rm**3)*(rum**4))
                jac_2(1,2,2) = (zum*((-2.*yum**4)-((xum**2)*(yum**2))
     1                         +(xum**4)+((xum**2)*(zum**2))))
     2                         /((cy_rm**3)*(rum**4))
                jac_2(1,2,3) = 2.*zum*(cy_rm)/(rum**4)
                jac_2(1,3,1) = (2*xum*yum)/((xum**2+yum**2)**2)
                jac_2(1,3,2) = -(2*xum*yum)/((xum**2+yum**2)**2)
                jac_2(1,3,3) = 0.

                jac_2(2,1,1) = (yua**2+zua**2)/(rua**3)
                jac_2(2,1,2) = (xua**2+zua**2)/(rua**3)
                jac_2(2,1,3) = (xua**2+yua**2)/(rua**3)
                jac_2(2,2,1) = (zua*((-2.*xua**4)-((xua**2)*(yua**2))
     1                         +(yua**4)+((yua**2)*(zua**2))))
     2                         /((cy_ra**3)*(rua**4))
                jac_2(2,2,2) = (zua*((-2.*yua**4)-((xua**2)*(yua**2))
     1                         +(xua**4)+((xua**2)*(zua**2))))
     2                         /((cy_ra**3)*(rua**4))
                jac_2(1,2,3) = 2.*zua*(cy_ra)/(rua**4)
                jac_2(2,3,1) = (2*xua*yua)/((xua**2+yua**2)**2)
                jac_2(2,3,2) = -(2*xua*yua)/((xua**2+yua**2)**2)
                jac_2(2,3,3) = 0.
! Spherical coordinates in terms of non-boosted
! cartesian coordinates                

                theta_m = acos(zum/rum)
                theta_a = acos(zua/rua)
! Since y and x coordinates are the same for the monpole and
! antimonopole, we use monopole-centered coordinates here to
! define the azimuthal angle.

                sph_phi = atan2(yum,xum)

! dimensionless distances (in units of vector mass):
!Need to correct for su2u1!

                sm=vev*rum
                sa=vev*rua

! -- end: setup variables ------


       S2=sin(theta_m/2.)
       S2b=sin(theta_a/2.)
       C2=cos(theta_m/2.)
       C2b=cos(theta_a/2.)

       dxhatn(1,i,j,k)=0.5*(((cmplx(cos(twist),sin(twist)))
     1     *(C2*S2b*(jac(1,2,1))+S2*C2b*(jac(2,2,1))))
     2     -(S2*C2b*(jac(1,2,1)))-(C2*S2b*(jac(2,2,1))))

       dxhatn(2,i,j,k)=0.5*((cmplx(cos(sph_phi),sin(sph_phi)))*(
     1     (C2*C2b*(jac(1,2,1)))-(S2*S2b*(jac(2,2,1)))
     2     +(2.*S2*C2b*(cmplx(0.,1.))*(jac(1,3,1)))
     3     +((cmplx(cos(twist),-sin(twist)))
     4     *((S2*S2b*(jac(1,2,1)))-(C2*C2b*(jac(2,2,1)))
     5     -(2.*C2*S2b*(cmplx(0.,1.))*(jac(1,3,1)))))))

       dyhatn(1,i,j,k)=0.5*(((cmplx(cos(twist),sin(twist)))
     1     *(C2*S2b*(jac(1,2,2))+S2*C2b*(jac(2,2,2))))
     2     -(S2*C2b*(jac(1,2,2)))-(C2*S2b*(jac(2,2,2))))

       dyhatn(2,i,j,k)=0.5*((cmplx(cos(sph_phi),sin(sph_phi)))*(
     1     (C2*C2b*(jac(1,2,2)))-(S2*S2b*(jac(2,2,2)))
     2     +(2.*S2*C2b*(cmplx(0.,1.))*(jac(1,3,2)))
     3     +((cmplx(cos(twist),-sin(twist)))
     4     *((S2*S2b*(jac(1,2,2)))-(C2*C2b*(jac(2,2,2)))
     5     -(2.*C2*S2b*(cmplx(0.,1.))*(jac(1,3,2)))))))

       dzhatn(1,i,j,k)=0.5*(((cmplx(cos(twist),sin(twist)))
     1     *(C2*S2b*(jac(1,2,3))+S2*C2b*(jac(2,2,3))))
     2     -(S2*C2b*(jac(1,2,3)))-(C2*S2b*(jac(2,2,3))))

       dzhatn(2,i,j,k)=0.5*((cmplx(cos(sph_phi),sin(sph_phi)))*(
     1     (C2*C2b*(jac(1,2,3)))-(S2*S2b*(jac(2,2,3)))
     2     +(2.*S2*C2b*(cmplx(0.,1.))*(jac(1,3,3)))
     3     +((cmplx(cos(twist),-sin(twist)))
     4     *((S2*S2b*(jac(1,2,3)))-(C2*C2b*(jac(2,2,3)))
     5     -(2.*C2*S2b*(cmplx(0.,1.))*(jac(1,3,3)))))))

       d2xhatn(1,i,j,k)=0.25*(((cmplx(cos(twist),sin(twist)))
     1    *(-(S2*S2b*((jac(1,2,1))**2))-(S2*S2b*((jac(2,2,1))**2))
     2     +(2.*(C2*C2b*(jac(1,2,1))*jac(2,2,1)))
     3     +(2.*C2*S2b*(jac_2(1,2,1)))+(2.*S2*C2b*(jac_2(2,2,1)))))
     4     -(C2*C2b*((jac(1,2,1))**2))-(C2*C2b*((jac(2,2,1))**2))
     5     +(2.*S2*S2b*(jac(1,2,1))*(jac(2,2,1)))
     6     -(S2*C2b*(jac_2(1,2,1)))-(C2*S2b*(jac_2(2,2,1))))

       d2xhatn(2,i,j,k)=((cmplx(0.,1.))*(dxhatn(2,i,j,k))*(jac(1,3,1)))
     1     +(0.25*(cmplx(cos(sph_phi),sin(sph_phi)))*(
     2     -(S2*C2b*(((jac(1,2,1))**2)+((jac(2,2,1))**2)))
     3     -(2.*C2*S2b*(jac(1,2,1))*(jac(2,2,1)))
     4     +(2.*C2*C2b*(jac_2(1,2,1)))-(2.*S2*S2b*(jac_2(2,2,1)))
     5     +((cmplx(0.,1.))*(((jac(1,3,1))*((2.*C2*C2b*(jac(1,2,1)))
     6     -(2.*S2*C2b*(jac(2,2,1)))))+(4.*S2*C2b*(jac(1,3,1)))))
     7     +((cmplx(cos(twist),-sin(twist)))*(
     8     (C2*S2b*(((jac(1,2,1))**2)+((jac(2,2,1))**2)))
     9     +(2.*S2*C2b*(jac(1,2,1))*(jac(2,2,1)))
     1     +(2.*S2*S2b*(jac_2(1,2,1)))-(2.*C2*C2b*(jac_2(2,2,1)))
     2     +((cmplx(0.,1.))*(((jac(1,3,1))*((2.*S2*S2b*(jac(1,2,1)))
     3     -(2.*C2*C2b*(jac(2,2,1)))))-(4.*C2*S2b*(jac(1,3,1)))))))))
          

       d2yhatn(1,i,j,k)=0.25*(((cmplx(cos(twist),sin(twist)))
     1    *(-(S2*S2b*((jac(1,2,2))**2))-(S2*S2b*((jac(2,2,2))**2))
     2     +(2.*(C2*C2b*(jac(1,2,2))*jac(2,2,2)))
     3     +(2.*C2*S2b*(jac_2(1,2,2)))+(2.*S2*C2b*(jac_2(2,2,2)))))
     4     -(C2*C2b*((jac(1,2,2))**2))-(C2*C2b*((jac(2,2,2))**2))
     5     +(2.*S2*S2b*(jac(1,2,2))*(jac(2,2,2)))
     6     -(S2*C2b*(jac_2(1,2,2)))-(C2*S2b*(jac_2(2,2,2))))

       d2yhatn(2,i,j,k)=((cmplx(0.,1.))*(dyhatn(2,i,j,k))*(jac(1,3,2)))
     1     +(0.25*(cmplx(cos(sph_phi),sin(sph_phi)))*(
     2     -(S2*C2b*(((jac(1,2,2))**2)+((jac(2,2,2))**2)))
     3     -(2.*C2*S2b*(jac(1,2,2))*(jac(2,2,2)))
     4     +(2.*C2*C2b*(jac_2(1,2,2)))-(2.*S2*S2b*(jac_2(2,2,2)))
     5     +((cmplx(0.,1.))*(((jac(1,3,2))*((2.*C2*C2b*(jac(1,2,2)))
     6     -(2.*S2*C2b*(jac(2,2,2)))))+(4.*S2*C2b*(jac(1,3,2)))))
     7     +((cmplx(cos(twist),-sin(twist)))*(
     8     (C2*S2b*(((jac(1,2,2))**2)+((jac(2,2,2))**2)))
     9     +(2.*S2*C2b*(jac(1,2,2))*(jac(2,2,2)))
     1     +(2.*S2*S2b*(jac_2(1,2,2)))-(2.*C2*C2b*(jac_2(2,2,2)))
     2     +((cmplx(0.,1.))*(((jac(1,3,2))*((2.*S2*S2b*(jac(1,2,2)))
     3     -(2.*C2*C2b*(jac(2,2,2)))))-(4.*C2*S2b*(jac(1,3,2)))))))))

       d2zhatn(1,i,j,k)=0.25*(((cmplx(cos(twist),sin(twist)))
     1    *(-(S2*S2b*((jac(1,2,3))**2))-(S2*S2b*((jac(2,2,3))**2))
     2     +(2.*(C2*C2b*(jac(1,2,3))*jac(2,2,3)))
     3     +(2.*C2*S2b*(jac_2(1,2,3)))+(2.*S2*C2b*(jac_2(2,2,3)))))
     4     -(C2*C2b*((jac(1,2,3))**2))-(C2*C2b*((jac(2,2,3))**2))
     5     +(2.*S2*S2b*(jac(1,2,3))*(jac(2,2,3)))
     6     -(S2*C2b*(jac_2(1,2,3)))-(C2*S2b*(jac_2(2,2,3))))

       d2zhatn(2,i,j,k)=((cmplx(0.,1.))*(dzhatn(2,i,j,k))*(jac(1,3,3)))
     1     +(0.25*(cmplx(cos(sph_phi),sin(sph_phi)))*(
     2     -(S2*C2b*(((jac(1,2,3))**2)+((jac(2,2,3))**2)))
     3     -(2.*C2*S2b*(jac(1,2,3))*(jac(2,2,3)))
     4     +(2.*C2*C2b*(jac_2(1,2,3)))-(2.*S2*S2b*(jac_2(2,2,3)))
     5     +((cmplx(0.,1.))*(((jac(1,3,3))*((2.*C2*C2b*(jac(1,2,3)))
     6     -(2.*S2*C2b*(jac(2,2,3)))))+(4.*S2*C2b*(jac(1,3,3)))))
     7     +((cmplx(cos(twist),-sin(twist)))*(
     8     (C2*S2b*(((jac(1,2,3))**2)+((jac(2,2,3))**2)))
     9     +(2.*S2*C2b*(jac(1,2,3))*(jac(2,2,3)))
     1     +(2.*S2*S2b*(jac_2(1,2,3)))-(2.*C2*C2b*(jac_2(2,2,3)))
     2     +((cmplx(0.,1.))*(((jac(1,3,3))*((2.*S2*S2b*(jac(1,2,3)))
     3     -(2.*C2*C2b*(jac(2,2,3)))))-(4.*C2*S2b*(jac(1,3,3)))))))))

                 enddo
             enddo
         enddo


!-----------end derivatives hatn--------!



!=======================================================================
! now for the initial gauge fields and the electric fields:
!=======================================================================

        do k=1,localsize
            do j=1,localsize
                do i=1,localsize
 
! -- begin: setup variables ------
! (x,y,z) for lattice site:
  
                iglobal = i+lxb-1
                jglobal = j+lyl-1
                kglobal = k+lzd-1

                x=float(iglobal)*dx
                y=float(jglobal)*dx
                z=float(kglobal)*dx

! Unboosted coordinates

                xum = (x-xm)
                yum = (y-ym)
                zum = (z-zm)
                xua = (x-xa)
                yua = (y-ya)
                zua = (z-za)
                rum = sqrt(xum**2+yum**2+zum**2)
                rua = sqrt(xua**2+yua**2+zua**2)
                cy_rm = sqrt(xum**2+yum**2)
                cy_ra = sqrt(xua**2+yua**2)

! Spherical coordinates in terms of non-boosted
! cartesian coordinates                

                theta_m = acos(zum/rum)
                theta_a = acos(zua/rua)
! Since y and x coordinates are the same for the monpole and
! antimonopole, we use monopole-centered coordinates here to
! define the azimuthal angle.

                sph_phi = atan2(yum,xum)

! dimensionless distances (in units of vector mass):

                sm=vev*rum
                sa=vev*rua

! Jacobian elements
! First index:1/2 -> monopole/antimonopole centered coordinates
! These are the elements corresponding to the partial derivatives of
! the spherical coordinates w.r.t the un-boosted cartesian coordinates.

                jac(1,1,1) = xum/rum
                jac(1,1,2) = yum/rum
                jac(1,1,3) = zum/rum
                jac(1,2,1) = xum*zum/(rum**2*sqrt(xum**2+yum**2))
                jac(1,2,2) = yum*zum/(rum**2*sqrt(xum**2+yum**2))
                jac(1,2,3) = -sqrt(xum**2+yum**2)/(rum**2)
                jac(1,3,1) = -yum/(xum**2+yum**2)
                jac(1,3,2) = xum/(xum**2+yum**2)
                jac(1,3,3) = 0.

                jac(2,1,1) = xua/rua
                jac(2,1,2) = yua/rua
                jac(2,1,3) = zua/rua
                jac(2,2,1) = xua*zua/(rua**2*sqrt(xua**2+yua**2))
                jac(2,2,2) = yua*zua/(rua**2*sqrt(xua**2+yua**2))
                jac(2,2,3) = -sqrt(xua**2+yua**2)/(rua**2)
                jac(2,3,1) = -yua/(xua**2+yua**2)
                jac(2,3,2) = xua/(xua**2+yua**2)
                jac(2,3,3) = 0.

                
                jac_2(1,1,1) = (yum**2+zum**2)/(rum**3)
                jac_2(1,1,2) = (xum**2+zum**2)/(rum**3)
                jac_2(1,1,3) = (xum**2+yum**2)/(rum**3)
                jac_2(1,2,1) = (zum*((-2.*xum**4)-((xum**2)*(yum**2))
     1                         +(yum**4)+((yum**2)*(zum**2))))
     2                         /((cy_rm**3)*(rum**4))
                jac_2(1,2,2) = (zum*((-2.*yum**4)-((xum**2)*(yum**2))
     1                         +(xum**4)+((xum**2)*(zum**2))))
     2                         /((cy_rm**3)*(rum**4))
                jac_2(1,2,3) = 2.*zum*(cy_rm)/(rum**4)
                jac_2(1,3,1) = (2*xum*yum)/((xum**2+yum**2)**2)
                jac_2(1,3,2) = -(2*xum*yum)/((xum**2+yum**2)**2)
                jac_2(1,3,3) = 0.

                jac_2(2,1,1) = (yua**2+zua**2)/(rua**3)
                jac_2(2,1,2) = (xua**2+zua**2)/(rua**3)
                jac_2(2,1,3) = (xua**2+yua**2)/(rua**3)
                jac_2(2,2,1) = (zua*((-2.*xua**4)-((xua**2)*(yua**2))
     1                         +(yua**4)+((yua**2)*(zua**2))))
     2                         /((cy_ra**3)*(rua**4))
                jac_2(2,2,2) = (zua*((-2.*yua**4)-((xua**2)*(yua**2))
     1                         +(xua**4)+((xua**2)*(zua**2))))
     2                         /((cy_ra**3)*(rua**4))
                jac_2(2,2,3) = 2.*zua*(cy_ra)/(rua**4)
                jac_2(2,3,1) = (2*xua*yua)/((xua**2+yua**2)**2)
                jac_2(2,3,2) = -(2*xua*yua)/((xua**2+yua**2)**2)
                jac_2(2,3,3) = 0.
                
!-----------------------------------------------------------------------
! monopole and antimonopole gauge profile function (1-K(r)):
! Note -- no factor of 1/r as in icMMbar.f and icMMbarTwisted.f
! Correction factors for profiles as evaluated by Ayush Saurabh:
!-----------------------------------------------------------------------

        correctionm=(1.+(1.-sqrt(lambda))*sm**2/4.+sm**4/16.)/
     1                (1.+sm**2/4.+sm**4/16.)
        correctiona=(1.+(1.-sqrt(lambda))*sa**2/4.+sa**4/16.)/
     1                (1.+sa**2/4.+sa**4/16.)
        if(sm.ne.0.) then
        wprofilem=(1./tanh(sm)-(1.+ms*sm)*exp(-ms*sm)/sm)/gw
        yprofilem=(1./tanh(sm)-(1.+ms*sm)*exp(-ms*sm)/sm)/gy
        else
        wprofilem=0.
        yprofilem=0.
        endif
        
        if(sa.ne.0.) then
        wprofilea=(1./tanh(sa)-(1.+ms*sa)*exp(-ms*sa)/sa)/gw
        yprofilea=(1./tanh(sa)-(1.+ms*sa)*exp(-ms*sa)/sa)/gy
        else
        wprofilea=0.
        yprofilea=0.
        endif

!        sprofile=(-(((tanh(-z+zm+szoff)+1.)/2.
!     1            +(tanh(z-za+szoff)+1.)/2.-1.)/(tanh(zm+szoff)))
!     2            *(cy_ra/sinh(cy_ra))+1.)

        sprofile=(-1.*((0.5*(tanh(-z+zm+szoff))
     1                  +0.5*(tanh(z-za+szoff)))/(tanh(zm+szoff)))
     2                  *(exp(-(cy_ra**2)*vev**2))+1.)

! square of scalar field su2u1:

       phisq=f(1,i,j,k)**2+f(2,i,j,k)**2+f(3,i,j,k)**2+f(4,i,j,k)**2

       if(phisq.ne.0.) then

! Define the double sided derivative term $\phi^\dag(\partial_\mu\phi)
! -(\partial_\mu\phi)^\dag\phi)$
!      if((cy_ra*cy_rm).ne.0.) then

       db_phi_x = conjg(hatn(1,i,j,k))*dxhatn(1,i,j,k)
     1           +conjg(hatn(2,i,j,k))*dxhatn(2,i,j,k)
     2           -conjg(dxhatn(1,i,j,k))*hatn(1,i,j,k)
     3           -conjg(dxhatn(2,i,j,k))*hatn(2,i,j,k)
       
       db_phi_y = conjg(hatn(1,i,j,k))*dyhatn(1,i,j,k)
     1           +conjg(hatn(2,i,j,k))*dyhatn(2,i,j,k)
     2           -conjg(dyhatn(1,i,j,k))*hatn(1,i,j,k)
     3           -conjg(dyhatn(2,i,j,k))*hatn(2,i,j,k)

       db_phi_z = conjg(hatn(1,i,j,k))*dzhatn(1,i,j,k)
     1           +conjg(hatn(2,i,j,k))*dzhatn(2,i,j,k)
     2           -conjg(dzhatn(1,i,j,k))*hatn(1,i,j,k)
     3           -conjg(dzhatn(2,i,j,k))*hatn(2,i,j,k)

!---Omitting the following expressions at \rho_m or \rho_a=0
!---Because of singularities

!      else
!        db_phi_x=cmplx(0.,0.)
!        db_phi_y=cmplx(0.,0.)
!        db_phi_z=cmplx(0.,0.)
!      endif

!      if((cy_rm*cy_ra*rua*rum).ne.0.)then
       n_1=(cos(twist-sph_phi)*(
     1     -ctw*cos(theta_a)*sin(theta_m)+cos(theta_m)*sin(theta_a))
     2     -stw*sin(theta_m)*sin(twist-sph_phi))
       
       n_2=(-cos(twist-sph_phi)*stw*sin(theta_m)+(
     1      ctw*cos(theta_a)*sin(theta_m)-cos(theta_m)*sin(theta_a))
     2      *sin(twist-sph_phi))

       n_3=(-cos(theta_m)*cos(theta_a)-ctw*sin(theta_m)*sin(theta_a))

       n_1_x=((-jac(1,2,1)*(cos(twist-sph_phi)*(
     1    (ctw*cos(theta_m)*cos(theta_a))+(sin(theta_m)*sin(theta_a)))
     2   +(cos(theta_m)*stw*sin(twist-sph_phi))))
     3   +((cos(twist-sph_phi)*jac(2,2,1))*(
     4    (cos(theta_m)*cos(theta_a))+(ctw*sin(theta_m)*sin(theta_a))))
     5   +(jac(1,3,1)*(
     6    (stw*cos(twist-sph_phi)*sin(theta_m))+sin(twist-sph_phi)*(
     7 -(ctw*cos(theta_a)*sin(theta_m))+(cos(theta_m)*sin(theta_a))))))

       n_1_y=((-jac(1,2,2)*(cos(twist-sph_phi)*(
     1    (ctw*cos(theta_m)*cos(theta_a))+(sin(theta_m)*sin(theta_a)))
     2   +(cos(theta_m)*stw*sin(twist-sph_phi))))
     3   +((cos(twist-sph_phi)*jac(2,2,2))*(
     4    (cos(theta_m)*cos(theta_a))+(ctw*sin(theta_m)*sin(theta_a))))
     5   +(jac(1,3,2)*(
     6    (stw*cos(twist-sph_phi)*sin(theta_m))+sin(twist-sph_phi)*(
     7 -(ctw*cos(theta_a)*sin(theta_m))+(cos(theta_m)*sin(theta_a))))))

       n_1_z=((-jac(1,2,3)*(cos(twist-sph_phi)*(
     1    (ctw*cos(theta_m)*cos(theta_a))+(sin(theta_m)*sin(theta_a)))
     2   +(cos(theta_m)*stw*sin(twist-sph_phi))))
     3   +((cos(twist-sph_phi)*jac(2,2,3))*(
     4    (cos(theta_m)*cos(theta_a))+(ctw*sin(theta_m)*sin(theta_a))))
     5   +(jac(1,3,3)*(
     6    (stw*cos(twist-sph_phi)*sin(theta_m))+sin(twist-sph_phi)*(
     7 -(ctw*cos(theta_a)*sin(theta_m))+(cos(theta_m)*sin(theta_a))))))

       n_2_x=((jac(1,2,1)*((-cos(theta_m)*stw*cos(twist-sph_phi))
     1  +(sin(twist-sph_phi)*(
     2  (ctw*cos(theta_m)*cos(theta_a))+(sin(theta_m)*sin(theta_a))))))
     3  -((sin(twist-sph_phi)*jac(2,2,1))*(
     4  (cos(theta_m)*cos(theta_a))+(ctw*sin(theta_m)*sin(theta_a))))
     5  -(jac(1,3,1)*(cos(twist-sph_phi)*(
     6  (ctw*cos(theta_a)*sin(theta_m))-(cos(theta_m)*sin(theta_a)))
     7  +(stw*sin(theta_m)*sin(twist-sph_phi)))))

       n_2_y=((jac(1,2,2)*((-cos(theta_m)*stw*cos(twist-sph_phi))
     1  +(sin(twist-sph_phi)*(
     2  (ctw*cos(theta_m)*cos(theta_a))+(sin(theta_m)*sin(theta_a))))))
     3  -((sin(twist-sph_phi)*jac(2,2,2))*(
     4  (cos(theta_m)*cos(theta_a))+(ctw*sin(theta_m)*sin(theta_a))))
     5  -(jac(1,3,2)*(cos(twist-sph_phi)*(
     6  (ctw*cos(theta_a)*sin(theta_m))-(cos(theta_m)*sin(theta_a)))
     7  +(stw*sin(theta_m)*sin(twist-sph_phi)))))

       n_2_z=((jac(1,2,3)*((-cos(theta_m)*stw*cos(twist-sph_phi))
     1  +(sin(twist-sph_phi)*(
     2  (ctw*cos(theta_m)*cos(theta_a))+(sin(theta_m)*sin(theta_a))))))
     3  -((sin(twist-sph_phi)*jac(2,2,3))*(
     4  (cos(theta_m)*cos(theta_a))+(ctw*sin(theta_m)*sin(theta_a))))
     5  -(jac(1,3,3)*(cos(twist-sph_phi)*(
     6  (ctw*cos(theta_a)*sin(theta_m))-(cos(theta_m)*sin(theta_a)))
     7  +(stw*sin(theta_m)*sin(twist-sph_phi)))))

       n_3_x=((jac(1,2,1)*(
     1   (cos(theta_a)*sin(theta_m))-(ctw*cos(theta_m)*sin(theta_a))))
     2  +(jac(2,2,1)*(
     3  -(ctw*cos(theta_a)*sin(theta_m))+(cos(theta_m)*sin(theta_a)))))

       n_3_y=((jac(1,2,2)*(
     1   (cos(theta_a)*sin(theta_m))-(ctw*cos(theta_m)*sin(theta_a))))
     2  +(jac(2,2,2)*(
     3  -(ctw*cos(theta_a)*sin(theta_m))+(cos(theta_m)*sin(theta_a)))))

       n_3_z=((jac(1,2,3)*(
     1   (cos(theta_a)*sin(theta_m))-(ctw*cos(theta_m)*sin(theta_a))))
     2  +(jac(2,2,3)*(
     3  -(ctw*cos(theta_a)*sin(theta_m))+(cos(theta_m)*sin(theta_a)))))


! W_\mu^1:
       f(5,i,j,k)=0.

       f(6,i,j,k)=gw*sprofile*wprofilem*wprofilea*(
     1            ((REAl(cmplx(0.,1.)*n_1*db_phi_x)*(cw**2))/(vev**2))
     2           -n_2*n_3_x+n_3*n_2_x)

       f(7,i,j,k)=gw*sprofile*wprofilem*wprofilea*(
     1            ((REAL(cmplx(0.,1.)*n_1*db_phi_y)*(cw**2))/(vev**2))
     2           -n_2*n_3_y+n_3*n_2_y)

       f(8,i,j,k)=gw*sprofile*wprofilem*wprofilea*(
     1            ((REAL(cmplx(0.,1.)*n_1*db_phi_z)*(cw**2))/(vev**2))
     2           -n_2*n_3_z+n_3*n_2_z)

! W_\mu^2:
       f(9,i,j,k)=0.

       f(10,i,j,k)=gw*sprofile*wprofilem*wprofilea*(
     1            ((REAL(cmplx(0.,1.)*n_2*db_phi_x)*(cw**2))/(vev**2))
     2           -n_3*n_1_x+n_1*n_3_x)

       f(11,i,j,k)=gw*sprofile*wprofilem*wprofilea*(
     1            ((REAL(cmplx(0.,1.)*n_2*db_phi_y)*(cw**2))/(vev**2))
     2           -n_3*n_1_y+n_1*n_3_y)

       f(12,i,j,k)=gw*sprofile*wprofilem*wprofilea*(
     1            ((REAL(cmplx(0.,1.)*n_2*db_phi_z)*(cw**2))/(vev**2))
     2           -n_3*n_1_z+n_1*n_3_z)

! W_\mu^3:
       f(13,i,j,k)=0.

       f(14,i,j,k)=gw*sprofile*wprofilem*wprofilea*(
     1            ((REAL(cmplx(0.,1.)*n_3*db_phi_x)*(cw**2))/(vev**2))
     2           -n_1*n_2_x+n_2*n_1_x)

       f(15,i,j,k)=gw*sprofile*wprofilem*wprofilea*(
     1            ((REAL(cmplx(0.,1.)*n_3*db_phi_y)*(cw**2))/(vev**2))
     2           -n_1*n_2_y+n_2*n_1_y)

       f(16,i,j,k)=gw*sprofile*wprofilem*wprofilea*(
     1            ((REAL(cmplx(0.,1.)*n_3*db_phi_z)*(cw**2))/(vev**2))
     2           -n_1*n_2_z+n_2*n_1_z)

! B_0 = 0:
       f(17,i,j,k)=0.

! B_1:
       f(18,i,j,k)=-gy*yprofilem*yprofilea*REAL(cmplx(0.,1.)*db_phi_x)
     1             *(sw**2)*sprofile/(vev**2)
! B_2
       f(19,i,j,k)=-gy*yprofilem*yprofilea*REAL(cmplx(0.,1.)*db_phi_y)
     1             *(sw**2)*sprofile/(vev**2)
! B_3
       f(20,i,j,k)=-gy*yprofilem*yprofilea*REAL(cmplx(0.,1.)*db_phi_z)
     1             *(sw**2)*sprofile/(vev**2)

       else
       f(5,i,j,k)=0.
       f(6,i,j,k)=0.
       f(7,i,j,k)=0.
       f(8,i,j,k)=0.
       f(9,i,j,k)=0.
       f(10,i,j,k)=0.
       f(11,i,j,k)=0.
       f(12,i,j,k)=0.
       f(13,i,j,k)=0.
       f(14,i,j,k)=0.
       f(15,i,j,k)=0.
       f(16,i,j,k)=0.
       f(17,i,j,k)=0.
       f(18,i,j,k)=0.
       f(19,i,j,k)=0.
       f(20,i,j,k)=0.

       endif

                enddo
            enddo
        enddo

!=======================================================================
! gauge functions $\Gamma^a = \partial_i W^a_i$,$\Xi_i = \partial_iB_i$:
! This uses a call to derivatives. That is why it is necessary to
! evaluate the gauge functions in a separate do loop (after all the
! gauge fields have been specified).
!=======================================================================

!-------------Set local indices for gamma -x----------------------------

                bx=4
                fx=localsize-3

!-------------Set local indices for gamma -y----------------------------

                ly=4
                ry=localsize-3

!-------------Set local indices for gamma -z----------------------------

                dz=4
                uz=localsize-3

        call icgamma(f,bx,fx,ly,ry,dz,uz,lxb,lyl,lzd)


        return
        end
