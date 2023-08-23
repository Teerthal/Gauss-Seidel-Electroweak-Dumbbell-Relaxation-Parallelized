!=======================================================================
! Here we define our differencing scheme (some differences are higher
! order than others -- better to make everything the same order).
! Note: this routine is for relaxation codes. E.g. derivatives of 
! fd(n,i,j,k) ($\dot f$ in evolution codes) are not calculated.
!=======================================================================

!==============================================================================
!-6th order in lattice interior everywhere with 2nd Order in boundary layers--!
!==============================================================================

      subroutine derivatives(f,i,j,k,dfdx,dfdy,dfdz,
     1            d2fdx2,d2fdy2,d2fdz2,lxb,lyl,lzd)
      implicit none
      include 'parameters.inc'
      integer i,j,k,n,is,ilocal,jlocal,klocal,lxb,lyl,lzd
      real*8 f,fd
      real*8 dfdx,dfdy,dfdz,fs,d2fdx2,d2fdy2,d2fdz2

      dimension f(nf,1:localsize,1:localsize,1:localsize)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension d2fdx2(nf),d2fdy2(nf),d2fdz2(nf)

            ilocal = i+lxb-1
            jlocal = j+lyl-1
            klocal = k+lzd-1

      do n=1,nf

!=======================================================================
! x-derivatives:
!=======================================================================
! 6th order 1st- and 2nd- derivatives:

       if(abs(ilocal).le.latx-3) then

! central derivatives in lattice interior:

        dfdx(n)=(f(n,i+3,j,k)-9.*f(n,i+2,j,k)+45.*f(n,i+1,j,k)
     1   -45.*f(n,i-1,j,k)+9.*f(n,i-2,j,k)-f(n,i-3,j,k))/(60.*dx)
        d2fdx2(n)=(f(n,i+3,j,k)/90.-3.*f(n,i+2,j,k)/20.
     1   +3.*f(n,i+1,j,k)/2.-49.*f(n,i,j,k)/18.+3.*f(n,i-1,j,k)/2.
     2   -3.*f(n,i-2,j,k)/20.+f(n,i-3,j,k)/90.)/(dx**2)

       endif

! one-sided derivatives near boundaries:

        is=+1.
        if(ilocal.lt.0) is=-1

!c 6th order:
   !      if(abs(ilocal).eq.latx-2) then

   !      dfdx(n)=-float(is)*(f(n,i+is*2,j,k)/30.-2.*f(n,i+is*1,j,k)/5.
   !   1   -7.*f(n,i,j,k)/12.+4.*f(n,i-is*1,j,k)/3.-f(n,i-is*2,j,k)/2.
   !   2   +2.*f(n,i-is*3,j,k)/15.-f(n,i-is*4,j,k)/60.)/dx 
   !      d2fdx2(n)=(-13.*f(n,i+is*2,j,k)/180.+19.*f(n,i+is*1,j,k)/15.
   !   1   -7.*f(n,i,j,k)/3.+10.*f(n,i-is*1,j,k)/9.
   !   2   +f(n,i-is*2,j,k)/12.-f(n,i-is*3,j,k)/15.
   !   3   +f(n,i-is*4,j,k)/90.)/dx**2

   !      endif

!2nd order at boundary layer (latx-2)
        if(abs(ilocal).eq.latx-2) then
        dfdx(n)=-float(is)*(f(n,i+is*2,j,k)/12.-2.*f(n,i+is*1,j,k)/3.
     1 +0.*f(n,i,j,k)+2.*f(n,i-is*1,j,k)/3.-f(n,i-is*2,j,k)/12.)/dx

        d2fdx2(n)=(-1.*f(n,i+is*2,j,k)/12.+4.*f(n,i+is*1,j,k)/3.
     1   -5.*f(n,i,j,k)/2.+4.*f(n,i-is*1,j,k)/3.
     2   -f(n,i-is*2,j,k)/12.)/dx**2
        endif

!c 6th order:
   !      if(abs(ilocal).eq.latx-1) then

   !      dfdx(n)=-float(is)*(-f(n,i+is*1,j,k)/6.-77.*f(n,i,j,k)/60.
   !   1   +5.*f(n,i-is*1,j,k)/2.-5.*f(n,i-is*2,j,k)/3.
   !   2   +5.*f(n,i-is*3,j,k)/6.-f(n,i-is*4,j,k)/4.
   !   3   +f(n,i-is*5,j,k)/30.)/dx 
   !      d2fdx2(n)=(137.*f(n,i+is*1,j,k)/180.-49.*f(n,i,j,k)/60.
   !   1   -17.*f(n,i-is*1,j,k)/12.+47.*f(n,i-is*2,j,k)/18.
   !   2   -19.*f(n,i-is*3,j,k)/12.+31.*f(n,i-is*4,j,k)/60.
   !   3   -13.*f(n,i-is*5,j,k)/180.)/dx**2

      !   endif

!2nd order at boundary layer (latx-1)
      if(abs(ilocal).eq.latx-1) then

        dfdx(n)=-float(is)*(-f(n,i+is*1,j,k)/3.-f(n,i,j,k)/2.
     1   +f(n,i-is*1,j,k)-f(n,i-is*2,j,k)/6.)/dx 
        d2fdx2(n)=(f(n,i+is*1,j,k)-2.*f(n,i,j,k)
     1   +f(n,i-is*1,j,k)+0.*f(n,i-is*2,j,k))/dx**2

        endif
        
!        if(abs(i).eq.latx) then
!c 6th order:
!        dfdx(n)=-float(is)*(-49.*f(n,i,j,k)/20.+6.*f(n,i-is*1,j,k)
!     1   -15.*f(n,i-is*2,j,k)/2.+20.*f(n,i-is*3,j,k)/3.
!     2   -15.*f(n,i-is*4,j,k)/4.+6.*f(n,i-is*5,j,k)/5.
!     3   -f(n,i-is*6,j,k)/6.)/dx
!        dfddx(n)=-float(is)*(-49.*fd(n,i,j,k)/20.+6.*fd(n,i-is*1,j,k)
!     1   -15.*fd(n,i-is*2,j,k)/2.+20.*fd(n,i-is*3,j,k)/3.
!     2   -15.*fd(n,i-is*4,j,k)/4.+6.*fd(n,i-is*5,j,k)/5.
!     3   -fd(n,i-is*6,j,k)/6.)/dx
!        d2fdx2(n)=(203.*f(n,i,j,k)/45.-87.*f(n,i-is*1,j,k)/5.
!     1   +117.*f(n,i-is*2,j,k)/4.-254.*f(n,i-is*3,j,k)/9.
!     2   +33.*f(n,i-is*4,j,k)/2.-27.*f(n,i-is*5,j,k)/5.
!     3   +137.*f(n,i-is*6,j,k)/180.)/dx**2
!        endif

!         if(abs(ilocal).eq.latx) then

! ! 8th order:

!         dfdx(n)=-float(is)*(-761.*f(n,i,j,k)/280.+8.*f(n,i-is*1,j,k)
!      1   -14.*f(n,i-is*2,j,k)+56.*f(n,i-is*3,j,k)/3.
!      2   -35.*f(n,i-is*4,j,k)/2.+56.*f(n,i-is*5,j,k)/5.
!      3   -14.*f(n,i-is*6,j,k)/3.+8.*f(n,i-is*7,j,k)/7.
!      4   -f(n,i-is*8,j,k)/8.)/dx
!         d2fdx2(n)=(29531.*f(n,i,j,k)/5040.-962.*f(n,i-is*1,j,k)/35.
!      1   +621.*f(n,i-is*2,j,k)/10.-4006.*f(n,i-is*3,j,k)/45.
!      2   +691.*f(n,i-is*4,j,k)/8.-282.*f(n,i-is*5,j,k)/5.
!      3   +2143.*f(n,i-is*6,j,k)/90.-206.*f(n,i-is*7,j,k)/35.
!      4   +363.*f(n,i-is*8,j,k)/560.)/dx**2
!         endif

!        if(abs(ilocal).eq.latx) then
! c 4th order: (since one-sided higher order may be unstable -- see README)
!        dfdx(n)=-float(is)*(-25.*f(n,i,j,k)/12.+4.*f(n,i-is*1,j,k)
!      1   -3.*f(n,i-is*2,j,k)+4.*f(n,i-is*3,j,k)/3.
!      2   -f(n,i-is*4,j,k)/4.)/dx
!    !      dfddx(n)=-float(is)*(-25.*fd(n,i,j,k)/12.+4.*fd(n,i-is*1,j,k)
!    !   1   -3.*fd(n,i-is*2,j,k)+4.*fd(n,i-is*3,j,k)/3.
!    !   2   -fd(n,i-is*4,j,k)/4.)/dx
!         d2fdx2(n)=(35.*f(n,i,j,k)/12.-26.*f(n,i-is*1,j,k)/3.
!      1   +19.*f(n,i-is*2,j,k)/2.-14.*f(n,i-is*3,j,k)/3.
!      2   +11.*f(n,i-is*4,j,k)/12.)/dx**2
!         endif

!2nd order at boundary layer (latx)
       if(abs(ilocal).eq.latx) then
       dfdx(n)=-float(is)*(-3.*f(n,i,j,k)/2.+2.*f(n,i-is*1,j,k)
     1   -0.5*f(n,i-is*2,j,k))/dx

        d2fdx2(n)=(f(n,i,j,k)-2.*f(n,i-is*1,j,k)
     1   +f(n,i-is*2,j,k))/dx**2
        endif


!=======================================================================
! y-derivatives:
!=======================================================================
! 6th order 1st- and 2nd- derivatives:

       if(abs(jlocal).le.laty-3) then
! central derivatives in lattice interior:

        dfdy(n)=(f(n,i,j+3,k)-9.*f(n,i,j+2,k)+45.*f(n,i,j+1,k)
     1   -45.*f(n,i,j-1,k)+9.*f(n,i,j-2,k)-f(n,i,j-3,k))/(60.*dx)
        d2fdy2(n)=(f(n,i,j+3,k)/90.-3.*f(n,i,j+2,k)/20.
     1   +3.*f(n,i,j+1,k)/2.-49.*f(n,i,j,k)/18.+3.*f(n,i,j-1,k)/2.
     2   -3.*f(n,i,j-2,k)/20.+f(n,i,j-3,k)/90.)/(dx**2)
       endif


! one-sided derivatives near boundaries:

        is=+1
        if(jlocal.lt.0) is=-1

!c 6th order:
   !      if(abs(jlocal).eq.laty-2) then
   !      dfdy(n)=-float(is)*(f(n,i,j+is*2,k)/30.-2.*f(n,i,j+is*1,k)/5.
   !   1   -7.*f(n,i,j,k)/12.+4.*f(n,i,j-is*1,k)/3.-f(n,i,j-is*2,k)/2.
   !   2   +2.*f(n,i,j-is*3,k)/15.-f(n,i,j-is*4,k)/60.)/dx 
   !      d2fdy2(n)=(-13.*f(n,i,j+is*2,k)/180.+19.*f(n,i,j+is*1,k)/15.
   !   1   -7.*f(n,i,j,k)/3.+10.*f(n,i,j-is*1,k)/9.
   !   2   +f(n,i,j-is*2,k)/12.-f(n,i,j-is*3,k)/15.
   !   3   +f(n,i,j-is*4,k)/90.)/dx**2
   !      endif

!2nd order at boundary layer (laty-2)
        if(abs(jlocal).eq.laty-2) then
        dfdy(n)=-float(is)*(f(n,i,j+is*2,k)/12.-2.*f(n,i,j+is*1,k)/3.
     1 +0.*f(n,i,j,k)+2.*f(n,i,j-is*1,k)/3.-f(n,i,j-is*2,k)/12.)/dx

        d2fdy2(n)=(-1.*f(n,i,j+is*2,k)/12.+4.*f(n,i,j+is*1,k)/3.
     1   -5.*f(n,i,j,k)/2.+4.*f(n,i,j-is*1,k)/3.
     2   -f(n,i,j-is*2,k)/12.)/dx**2
        endif

!c 6th order:
   !      if(abs(jlocal).eq.laty-1) then
   !      dfdy(n)=-float(is)*(-f(n,i,j+is*1,k)/6.-77.*f(n,i,j,k)/60.
   !   1   +5.*f(n,i,j-is*1,k)/2.-5.*f(n,i,j-is*2,k)/3.
   !   2   +5.*f(n,i,j-is*3,k)/6.-f(n,i,j-is*4,k)/4.
   !   3   +f(n,i,j-is*5,k)/30.)/dx 
   !      d2fdy2(n)=(137.*f(n,i,j+is*1,k)/180.-49.*f(n,i,j,k)/60.
   !   1   -17.*f(n,i,j-is*1,k)/12.+47.*f(n,i,j-is*2,k)/18.
   !   2   -19.*f(n,i,j-is*3,k)/12.+31.*f(n,i,j-is*4,k)/60.
   !   3   -13.*f(n,i,j-is*5,k)/180.)/dx**2
   !      endif

!2nd order at boundary layer (laty-1)
      if(abs(jlocal).eq.laty-1) then

        dfdy(n)=-float(is)*(-f(n,i,j+is*1,k)/3.-f(n,i,j,k)/2.
     1   +f(n,i,j-is*1,k)-f(n,i,j-is*2,k)/6.)/dx 
        d2fdy2(n)=(f(n,i,j+is*1,k)-2.*f(n,i,j,k)
     1   +f(n,i,j-is*1,k)+0.*f(n,i,j-is*2,k))/dx**2

        endif

      !   if(abs(jlocal).eq.laty) then

! 8th order:

   !      dfdy(n)=-float(is)*(-761.*f(n,i,j,k)/280.+8.*f(n,i,j-is*1,k)
   !   1   -14.*f(n,i,j-is*2,k)+56.*f(n,i,j-is*3,k)/3.
   !   2   -35.*f(n,i,j-is*4,k)/2.+56.*f(n,i,j-is*5,k)/5.
   !   3   -14.*f(n,i,j-is*6,k)/3.+8.*f(n,i,j-is*7,k)/7.
   !   4   -f(n,i,j-is*8,k)/8.)/dx
   !      d2fdy2(n)=(29531.*f(n,i,j,k)/5040.-962.*f(n,i,j-is*1,k)/35.
   !   1   +621.*f(n,i,j-is*2,k)/10.-4006.*f(n,i,j-is*3,k)/45.
   !   2   +691.*f(n,i,j-is*4,k)/8.-282.*f(n,i,j-is*5,k)/5.
   !   3   +2143.*f(n,i,j-is*6,k)/90.-206.*f(n,i,j-is*7,k)/35.
   !   4   +363.*f(n,i,j-is*8,k)/560.)/dx**2
   !      endif

!
!        if(abs(jlocal).eq.laty) then
! c 4th order: (since one-sided higher order may be unstable -- see README)
!          dfdy(n)=-float(is)*(-25.*f(n,i,j,k)/12.+4.*f(n,i,j-is*1,k)
!      1   -3.*f(n,i,j-is*2,k)+4.*f(n,i,j-is*3,k)/3.
!      2   -f(n,i,j-is*4,k)/4.)/dx
!    !      dfddy(n)=-float(is)*(-25.*fd(n,i,j,k)/12.+4.*fd(n,i,j-is*1,k)
!    !   1   -3.*fd(n,i,j-is*2,k)+4.*fd(n,i,j-is*3,k)/3.
!    !   2   -fd(n,i,j-is*4,k)/4.)/dx
!         d2fdy2(n)=(35.*f(n,i,j,k)/12.-26.*f(n,i,j-is*1,k)/3.
!      1   +19.*f(n,i,j-is*2,k)/2.-14.*f(n,i,j-is*3,k)/3.
!      2   +11.*f(n,i,j-is*4,k)/12.)/dx**2
!          endif

!2nd order at boundary layer (laty)
       if(abs(jlocal).eq.laty) then
       dfdy(n)=-float(is)*(-3.*f(n,i,j,k)/2.+2.*f(n,i,j-is*1,k)
     1   -0.5*f(n,i,j-is*2,k))/dx

        d2fdy2(n)=(f(n,i,j,k)-2.*f(n,i,j-is*1,k)
     1   +f(n,i,j-is*2,k))/dx**2
        endif

!=======================================================================
! z-derivatives:
!=======================================================================
! 6th order 1st- and 2nd- derivatives:

       if(abs(klocal).le.latz-3) then

! central derivatives in lattice interior:

        dfdz(n)=(f(n,i,j,k+3)-9.*f(n,i,j,k+2)+45.*f(n,i,j,k+1)
     1   -45.*f(n,i,j,k-1)+9.*f(n,i,j,k-2)-f(n,i,j,k-3))/(60.*dx)
        d2fdz2(n)=(f(n,i,j,k+3)/90.-3.*f(n,i,j,k+2)/20.
     1   +3.*f(n,i,j,k+1)/2.-49.*f(n,i,j,k)/18.+3.*f(n,i,j,k-1)/2.
     2   -3.*f(n,i,j,k-2)/20.+f(n,i,j,k-3)/90.)/(dx**2)
       endif

! one-sided derivatives near boundaries:

        is=+1
        if(klocal.lt.0) is=-1

!c 6th order:
   !      if(abs(klocal).eq.latz-2) then
   !      dfdz(n)=-float(is)*(f(n,i,j,k+is*2)/30.-2.*f(n,i,j,k+is*1)/5.
   !   1   -7.*f(n,i,j,k)/12.+4.*f(n,i,j,k-is*1)/3.-f(n,i,j,k-is*2)/2.
   !   2   +2.*f(n,i,j,k-is*3)/15.-f(n,i,j,k-is*4)/60.)/dx 
   !      d2fdz2(n)=(-13.*f(n,i,j,k+is*2)/180.+19.*f(n,i,j,k+is*1)/15.
   !   1   -7.*f(n,i,j,k)/3.+10.*f(n,i,j,k-is*1)/9.
   !   2   +f(n,i,j,k-is*2)/12.-f(n,i,j,k-is*3)/15.
   !   3   +f(n,i,j,k-is*4)/90.)/dx**2
   !      endif

!2nd order at boundary layer (latz-2)
        if(abs(klocal).eq.latz-2) then
        dfdz(n)=-float(is)*(f(n,i,j,k+is*2)/12.-2.*f(n,i,j,k+is*1)/3.
     1 +0.*f(n,i,j,k)+2.*f(n,i,j,k-is*1)/3.-f(n,i,j,k-is*2)/12.)/dx

        d2fdz2(n)=(-1.*f(n,i,j,k+is*2)/12.+4.*f(n,i,j,k+is*1)/3.
     1   -5.*f(n,i,j,k)/2.+4.*f(n,i,j,k-is*1)/3.
     2   -f(n,i,j,k-is*2)/12.)/dx**2
        endif

!c 6th order:        
   !      if(abs(klocal).eq.latz-1) then
   !      dfdz(n)=-float(is)*(-f(n,i,j,k+is*1)/6.-77.*f(n,i,j,k)/60.
   !   1   +5.*f(n,i,j,k-is*1)/2.-5.*f(n,i,j,k-is*2)/3.
   !   2   +5.*f(n,i,j,k-is*3)/6.-f(n,i,j,k-is*4)/4.
   !   3   +f(n,i,j,k-is*5)/30.)/dx 
   !      d2fdz2(n)=(137.*f(n,i,j,k+is*1)/180.-49.*f(n,i,j,k)/60.
   !   1   -17.*f(n,i,j,k-is*1)/12.+47.*f(n,i,j,k-is*2)/18.
   !   2   -19.*f(n,i,j,k-is*3)/12.+31.*f(n,i,j,k-is*4)/60.
   !   3   -13.*f(n,i,j,k-is*5)/180.)/dx**2
   !      endif

!2nd order at boundary layer (latz-1)        
      if(abs(klocal).eq.latz-1) then

        dfdz(n)=-float(is)*(-f(n,i,j,k+is*1)/3.-f(n,i,j,k)/2.
     1   +f(n,i,j,k-is*1)-f(n,i,j,k-is*2)/6.)/dx 
        d2fdz2(n)=(f(n,i,j,k+is*1)-2.*f(n,i,j,k)
     1   +f(n,i,j,k-is*1)+0.*f(n,i,j,k-is*2))/dx**2

        endif


!         if(abs(klocal).eq.latz) then

! c 8th order:

!         dfdz(n)=-float(is)*(-761.*f(n,i,j,k)/280.+8.*f(n,i,j,k-is*1)
!      1   -14.*f(n,i,j,k-is*2)+56.*f(n,i,j,k-is*3)/3.
!      2   -35.*f(n,i,j,k-is*4)/2.+56.*f(n,i,j,k-is*5)/5.
!      3   -14.*f(n,i,j,k-is*6)/3.+8.*f(n,i,j,k-is*7)/7.
!      4   -f(n,i,j,k-is*8)/8.)/dx
!         d2fdz2(n)=(29531.*f(n,i,j,k)/5040.-962.*f(n,i,j,k-is*1)/35.
!      1   +621.*f(n,i,j,k-is*2)/10.-4006.*f(n,i,j,k-is*3)/45.
!      2   +691.*f(n,i,j,k-is*4)/8.-282.*f(n,i,j,k-is*5)/5.
!      3   +2143.*f(n,i,j,k-is*6)/90.-206.*f(n,i,j,k-is*7)/35.
!      4   +363.*f(n,i,j,k-is*8)/560.)/dx**2
!         endif

!        if(abs(klocal).eq.latz) then
! c 4th order: (since one-sided higher order may be unstable -- see README)
!        dfdz(n)=-float(is)*(-25.*f(n,i,j,k)/12.+4.*f(n,i,j,k-is*1)
!      1   -3.*f(n,i,j,k-is*2)+4.*f(n,i,j,k-is*3)/3.
!      2   -f(n,i,j,k-is*4)/4.)/dx
!    !      dfddz(n)=-float(is)*(-25.*fd(n,i,j,k)/12.+4.*fd(n,i,j,k-is*1)
!    !   1   -3.*fd(n,i,j,k-is*2)+4.*fd(n,i,j,k-is*3)/3.
!    !   2   -fd(n,i,j,k-is*4)/4.)/dx
!         d2fdz2(n)=(35.*f(n,i,j,k)/12.-26.*f(n,i,j,k-is*1)/3.
!      1   +19.*f(n,i,j,k-is*2)/2.-14.*f(n,i,j,k-is*3)/3.
!      2   +11.*f(n,i,j,k-is*4)/12.)/dx**2
!         endif

!2nd order at boundary layer (latz)
       if(abs(klocal).eq.latz) then
       dfdz(n)=-float(is)*(-3.*f(n,i,j,k)/2.+2.*f(n,i,j,k-is*1)
     1   -0.5*f(n,i,j,k-is*2))/dx

        d2fdz2(n)=(f(n,i,j,k)-2.*f(n,i,j,k-is*1)
     1   +f(n,i,j,k-is*2))/dx**2
        endif


       enddo

      return
      end
