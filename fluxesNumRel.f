!=========================================================
!The fluxes on the RHS of the Gauss Seidel relaxation
!called by euler evolve
!The equations for the nf=1:4 field laplacians and fluxes
!are stated but unnecessary since they are determined by 
!nf=25 component
!=========================================================
      subroutine fdflux(f,i,j,k,r,lxb,lyl,lzd,
     1                  hatn,dxhatn,dyhatn,dzhatn,
     2                  d2xhatn,d2yhatn,d2zhatn)

      implicit none
      include 'parameters.inc'
      integer i,j,k,n,nn,lxb,lyl,lzd,ilocal,jlocal,klocal
      real*8 f,r
      real*8 dfdx,dfdy,dfdz,fs,d2fdx2,d2fdy2,d2fdz2
      real*8 cd,g
      real*8 x,y,z,rr
      real*8 rmserror
      real*8 dfddx,dfddy,dfddz
      real*8 res
      real*8 laplacian
      real*8 im_tol

      complex*8 lap_mod_phi
      complex*8 hatn
      complex*8 dxhatn, dyhatn, dzhatn
      complex*8 d2xhatn, d2yhatn, d2zhatn
      complex*8 tem_mod_flx

      dimension f(nf,1:localsize,1:localsize,1:localsize)
      dimension r(nf)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension d2fdx2(nf),d2fdy2(nf),d2fdz2(nf)
      dimension laplacian(nf)
      dimension hatn(2,1:localsize,1:localsize,1:localsize)
      dimension dxhatn(2,1:localsize,1:localsize,1:localsize)
      dimension dyhatn(2,1:localsize,1:localsize,1:localsize)
      dimension dzhatn(2,1:localsize,1:localsize,1:localsize)
      dimension d2xhatn(2,1:localsize,1:localsize,1:localsize)
      dimension d2yhatn(2,1:localsize,1:localsize,1:localsize)
      dimension d2zhatn(2,1:localsize,1:localsize,1:localsize)
!=======================================================================
! fs=field strength; first index is the number of gauge fields.
!=======================================================================

      dimension fs(4,0:3,0:3),dfddx(nf),dfddy(nf),dfddz(nf)

!=======================================================================
! cd=covariant derivative; second index is the gauge field index.
!=======================================================================

      dimension cd(0:3,4),g(nf)
      dimension res(nf)

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

! Tolerance parameter for erroneous imaginary components of |\Phi|
! which should be strictly real

          im_tol=1e-10

!local here seem to imply global coordinates.
            ilocal = i+lxb-1
            jlocal = j+lyl-1
            klocal = k+lzd-1

! spatial derivatives to 4th order:

      call derivatives(f,i,j,k,dfdx,dfdy,dfdz,
     1                  d2fdx2,d2fdy2,d2fdz2,lxb,lyl,lzd)

! covariant derivatives:

      call covderiv(f,i,j,k,dfdx,dfdy,dfdz,cd)

! find gauge field strengths (uses derivatives):

      call fieldstrength(f,i,j,k,dfdx,dfdy,dfdz,fs)

! for convenience of writing:

      do n=1,nf
      g(n)=f(n,i,j,k)
      enddo

!=======================================================================
! define laplacian() such that the expression is accurate near the
! boundaries as well. For this relaxation code, we will use the boundary 
! conditions: $D_i\phi^a =0$ and $W^a_{ij}=0$. These relations then let
! us write second order derivatives in terms of first order derivatives.
! if on boundaries:
      
      if(abs(ilocal).eq.latx.or.abs(jlocal)
     1                          .eq.laty.or.abs(klocal).eq.latz) then

! laplacian on the boundary:
! scalar fields;
      laplacian(1)=-0.5*gw*((-g(6)*dfdx(4)-g(7)*dfdy(4)-g(8)*dfdz(4))
     2         -(-g(10)*dfdx(3)-g(11)*dfdy(3)-g(12)*dfdz(3))
     3         +(-g(14)*dfdx(2)-g(15)*dfdy(2)-g(16)*dfdz(2)))
     4 -0.5*gy*(-g(18)*dfdx(2)-g(19)*dfdy(2)-g(20)*dfdz(2))
     7 -0.5*((gw*g(23)+gy*g(24))*g(2)-gw*g(22)*g(3)+gw*g(21)*g(4))

      laplacian(2)=0.5*gw*((-g(6)*dfdx(3)-g(7)*dfdy(3)-g(8)*dfdz(3))
     2     +(-g(10)*dfdx(4)-g(11)*dfdy(4)-g(12)*dfdz(4))
     3     +(-g(14)*dfdx(1)-g(15)*dfdy(1)-g(16)*dfdz(1)))
     4 +0.5*gy*(-g(18)*dfdx(1)-g(19)*dfdy(1)-g(20)*dfdz(1))
     5 -0.5*((gw*g(23)+gy*g(24))*g(1)+gw*g(21)*g(3)+gw*g(22)*g(4))

      laplacian(3)=-0.5*gw*((-g(6)*dfdx(2)-g(7)*dfdy(2)-g(8)*dfdz(2))
     2     +(-g(10)*dfdx(1)-g(11)*dfdy(1)-g(12)*dfdz(1))
     3     -(-g(14)*dfdx(4)-g(15)*dfdy(4)-g(16)*dfdz(4)))
     4 -0.5*gy*(-g(18)*dfdx(4)-g(19)*dfdy(4)-g(20)*dfdz(4))
     5 +0.5*((-gw*g(23)+gy*g(24))*g(4)+gw*g(22)*g(1)+gw*g(21)*g(2))

      laplacian(4)=0.5*gw*((-g(6)*dfdx(1)-g(7)*dfdy(1)-g(8)*dfdz(1))
     2     -(-g(10)*dfdx(2)-g(11)*dfdy(2)-g(12)*dfdz(2))
     3     -(-g(14)*dfdx(3)-g(15)*dfdy(3)-g(16)*dfdz(3)))
     4 +0.5*gy*(-g(18)*dfdx(3)-g(19)*dfdy(3)-g(20)*dfdz(3))
     5 -0.5*((-gw*g(23)+gy*g(24))*g(3)+gw*g(21)*g(1)-gw*g(22)*g(2))

! gauge fields;
! W^1_i
          laplacian(5)=0.

      laplacian(6)=-(gw*( 
     2   -(dfdx(10)*f(14,i,j,k)-dfdx(14)*f(10,i,j,k))
     2   -(dfdy(10)*f(15,i,j,k)-dfdy(14)*f(11,i,j,k))
     2   -(dfdz(10)*f(16,i,j,k)-dfdz(14)*f(12,i,j,k)))
     5   -dfdx(21)-gw*(g(10)*g(23)-g(14)*g(22)))


      laplacian(7)=-(gw*( 
     2   -(dfdx(11)*f(14,i,j,k)-dfdx(15)*f(10,i,j,k))
     2   -(dfdy(11)*f(15,i,j,k)-dfdy(15)*f(11,i,j,k))
     2   -(dfdz(11)*f(16,i,j,k)-dfdz(15)*f(12,i,j,k)))
     5   -dfdy(21)-gw*(g(11)*g(23)-g(15)*g(22)))
          
      laplacian(8)=-(gw*( 
     2   -(dfdx(12)*f(14,i,j,k)-dfdx(16)*f(10,i,j,k))
     2   -(dfdy(12)*f(15,i,j,k)-dfdy(16)*f(11,i,j,k))
     2   -(dfdz(12)*f(16,i,j,k)-dfdz(16)*f(12,i,j,k)))
     5   -dfdz(21)-gw*(g(12)*g(23)-g(16)*g(22)))


! W^2_i
          laplacian(9)=0.

      laplacian(10)=-(gw*( 
     2   -(dfdx(14)*f(6,i,j,k)-dfdx(6)*f(14,i,j,k))
     2   -(dfdy(14)*f(7,i,j,k)-dfdy(6)*f(15,i,j,k))
     2   -(dfdz(14)*f(8,i,j,k)-dfdz(6)*f(16,i,j,k)))
     5   -dfdx(22)-gw*(g(14)*g(21)-g(6)*g(23)))

      laplacian(11)=-(gw*( 
     2   -(dfdx(15)*f(6,i,j,k)-dfdx(7)*f(14,i,j,k))
     2   -(dfdy(15)*f(7,i,j,k)-dfdy(7)*f(15,i,j,k))
     2   -(dfdz(15)*f(8,i,j,k)-dfdz(7)*f(16,i,j,k)))
     5   -dfdy(22)-gw*(g(15)*g(21)-g(7)*g(23)))

      laplacian(12)=-(gw*( 
     2   -(dfdx(16)*f(6,i,j,k)-dfdx(8)*f(14,i,j,k))
     2   -(dfdy(16)*f(7,i,j,k)-dfdy(8)*f(15,i,j,k))
     2   -(dfdz(16)*f(8,i,j,k)-dfdz(8)*f(16,i,j,k)))
     5   -dfdz(22)-gw*(g(16)*g(21)-g(8)*g(23)))

! W^3_i
          laplacian(13)=0.

      laplacian(14)=-(gw*( 
     2   -(dfdx(6)*f(10,i,j,k)-dfdx(10)*f(6,i,j,k))
     2   -(dfdy(6)*f(11,i,j,k)-dfdy(10)*f(7,i,j,k))
     2   -(dfdz(6)*f(12,i,j,k)-dfdz(10)*f(8,i,j,k)))
     5   -dfdx(23)-gw*(g(6)*g(22)-g(10)*g(21)))

      laplacian(15)=-(gw*( 
     2   -(dfdx(7)*f(10,i,j,k)-dfdx(11)*f(6,i,j,k))
     2   -(dfdy(7)*f(11,i,j,k)-dfdy(11)*f(7,i,j,k))
     2   -(dfdz(7)*f(12,i,j,k)-dfdz(11)*f(8,i,j,k)))
     5   -dfdy(23)-gw*(g(7)*g(22)-g(11)*g(21)))

      laplacian(16)=-(gw*( 
     2   -(dfdx(8)*f(10,i,j,k)-dfdx(12)*f(6,i,j,k))
     2   -(dfdy(8)*f(11,i,j,k)-dfdy(12)*f(7,i,j,k))
     2   -(dfdz(8)*f(12,i,j,k)-dfdz(12)*f(8,i,j,k)))
     5   -dfdz(23)-gw*(g(8)*g(22)-g(12)*g(21)))

!B_i
          laplacian(17)=0.

          laplacian(18)=dfdx(24)

          laplacian(19)=dfdy(24)

          laplacian(20)=dfdz(24)

! gauge constraints:

          laplacian(21)=0.
          laplacian(22)=0.
          laplacian(23)=0. 
          laplacian(24)=0.

          lap_mod_phi=(1./2.)*(
     1              conjg(d2xhatn(1,i,j,k))*cmplx(g(1),g(2))
     1             +conjg(d2xhatn(2,i,j,k))*cmplx(g(3),g(4))
     1             +conjg(d2yhatn(1,i,j,k))*cmplx(g(1),g(2))
     1             +conjg(d2yhatn(2,i,j,k))*cmplx(g(3),g(4))
     1             +conjg(d2zhatn(1,i,j,k))*cmplx(g(1),g(2))
     1             +conjg(d2zhatn(2,i,j,k))*cmplx(g(3),g(4))
     2             +conjg(cmplx(g(1),g(2)))*d2xhatn(1,i,j,k)
     2             +conjg(cmplx(g(3),g(4)))*d2xhatn(2,i,j,k)
     2             +conjg(cmplx(g(1),g(2)))*d2yhatn(1,i,j,k)
     2             +conjg(cmplx(g(3),g(4)))*d2yhatn(2,i,j,k)
     2             +conjg(cmplx(g(1),g(2)))*d2zhatn(1,i,j,k)
     2             +conjg(cmplx(g(3),g(4)))*d2zhatn(2,i,j,k)
     3             +conjg(hatn(1,i,j,k))
     3             *cmplx(laplacian(1),laplacian(2))
     3             +conjg(hatn(2,i,j,k))
     3             *cmplx(laplacian(3),laplacian(4))
     4             +conjg(cmplx(laplacian(1),laplacian(2)))
     4             *hatn(1,i,j,k)
     4             +conjg(cmplx(laplacian(3),laplacian(4)))
     4             *hatn(2,i,j,k)
     5             +2.*(conjg(dxhatn(1,i,j,k)))*cmplx(dfdx(1),dfdx(2))
     5             +2.*(conjg(dxhatn(2,i,j,k)))*cmplx(dfdx(3),dfdx(4))
     5             +2.*(conjg(dyhatn(1,i,j,k)))*cmplx(dfdy(1),dfdy(2))
     5             +2.*(conjg(dyhatn(2,i,j,k)))*cmplx(dfdy(3),dfdy(4))
     5             +2.*(conjg(dzhatn(1,i,j,k)))*cmplx(dfdz(1),dfdz(2))
     5             +2.*(conjg(dzhatn(2,i,j,k)))*cmplx(dfdz(3),dfdz(4))
     6             +2.*conjg(cmplx(dfdx(1),dfdx(2)))*dxhatn(1,i,j,k)
     6             +2.*conjg(cmplx(dfdx(3),dfdx(4)))*dxhatn(2,i,j,k)
     6             +2.*conjg(cmplx(dfdy(1),dfdy(2)))*dyhatn(1,i,j,k)
     6             +2.*conjg(cmplx(dfdy(3),dfdy(4)))*dyhatn(2,i,j,k)
     6             +2.*conjg(cmplx(dfdz(1),dfdz(2)))*dzhatn(1,i,j,k)
     6             +2.*conjg(cmplx(dfdz(3),dfdz(4)))*dzhatn(2,i,j,k))
          
! if not on boundaries:
      else
          do n=1,nf-1
          laplacian(n)=d2fdx2(n)+d2fdy2(n)+d2fdz2(n)
          enddo

          lap_mod_phi=d2fdx2(25)+d2fdy2(25)+d2fdz2(25)

      endif

!=======================================================================
! Scalar fluxes:
!=======================================================================

      r(1)=laplacian(1)
     1 -0.5*gw*((-g(6)*dfdx(4)-g(7)*dfdy(4)-g(8)*dfdz(4))
     2         -(-g(10)*dfdx(3)-g(11)*dfdy(3)-g(12)*dfdz(3))
     3         +(-g(14)*dfdx(2)-g(15)*dfdy(2)-g(16)*dfdz(2)))
     4 -0.5*gy*(-g(18)*dfdx(2)-g(19)*dfdy(2)-g(20)*dfdz(2))
     1 -0.5*gw*((-g(6)*cd(1,4)-g(7)*cd(2,4)-g(8)*cd(3,4))
     2     -(-g(10)*cd(1,3)-g(11)*cd(2,3)-g(12)*cd(3,3))
     3     +(-g(14)*cd(1,2)-g(15)*cd(2,2)-g(16)*cd(3,2)))
     4 -0.5*gy*(-g(18)*cd(1,2)-g(19)*cd(2,2)-g(20)*cd(3,2))
     1 -2.*lambda*(g(1)**2+g(2)**2+g(3)**2+g(4)**2-vev**2)*g(1)
     5 +0.5*((gw*g(23)+gy*g(24))*g(2)-gw*g(22)*g(3)+gw*g(21)*g(4))
c
      r(2)=laplacian(2)
     1 +0.5*gw*((-g(6)*dfdx(3)-g(7)*dfdy(3)-g(8)*dfdz(3))
     2     +(-g(10)*dfdx(4)-g(11)*dfdy(4)-g(12)*dfdz(4))
     3     +(-g(14)*dfdx(1)-g(15)*dfdy(1)-g(16)*dfdz(1)))
     4 +0.5*gy*(-g(18)*dfdx(1)-g(19)*dfdy(1)-g(20)*dfdz(1))
     1 +0.5*gw*((-g(6)*cd(1,3)-g(7)*cd(2,3)-g(8)*cd(3,3))
     2     +(-g(10)*cd(1,4)-g(11)*cd(2,4)-g(12)*cd(3,4))
     3     +(-g(14)*cd(1,1)-g(15)*cd(2,1)-g(16)*cd(3,1)))
     4 +0.5*gy*(-g(18)*cd(1,1)-g(19)*cd(2,1)-g(20)*cd(3,1))
     1 -2.*lambda*(g(1)**2+g(2)**2+g(3)**2+g(4)**2-vev**2)*g(2)
     5 -0.5*((gw*g(23)+gy*g(24))*g(1)+gw*g(21)*g(3)+gw*g(22)*g(4))
c
      r(3)=laplacian(3)
     1 -0.5*gw*((-g(6)*dfdx(2)-g(7)*dfdy(2)-g(8)*dfdz(2))
     2     +(-g(10)*dfdx(1)-g(11)*dfdy(1)-g(12)*dfdz(1))
     3     -(-g(14)*dfdx(4)-g(15)*dfdy(4)-g(16)*dfdz(4)))
     4 -0.5*gy*(-g(18)*dfdx(4)-g(19)*dfdy(4)-g(20)*dfdz(4))
     1 -0.5*gw*((-g(6)*cd(1,2)-g(7)*cd(2,2)-g(8)*cd(3,2))
     2     +(-g(10)*cd(1,1)-g(11)*cd(2,1)-g(12)*cd(3,1))
     3     -(-g(14)*cd(1,4)-g(15)*cd(2,4)-g(16)*cd(3,4)))
     4 -0.5*gy*(-g(18)*cd(1,4)-g(19)*cd(2,4)-g(20)*cd(3,4))
     1 -2.*lambda*(g(1)**2+g(2)**2+g(3)**2+g(4)**2-vev**2)*g(3)
     5 +0.5*((-gw*g(23)+gy*g(24))*g(4)+gw*g(22)*g(1)+gw*g(21)*g(2))
c
      r(4)=laplacian(4)
     1 +0.5*gw*((-g(6)*dfdx(1)-g(7)*dfdy(1)-g(8)*dfdz(1))
     2     -(-g(10)*dfdx(2)-g(11)*dfdy(2)-g(12)*dfdz(2))
     3     -(-g(14)*dfdx(3)-g(15)*dfdy(3)-g(16)*dfdz(3)))
     4 +0.5*gy*(-g(18)*dfdx(3)-g(19)*dfdy(3)-g(20)*dfdz(3))
     1 +0.5*gw*((-g(6)*cd(1,1)-g(7)*cd(2,1)-g(8)*cd(3,1))
     2     -(-g(10)*cd(1,2)-g(11)*cd(2,2)-g(12)*cd(3,2))
     3     -(-g(14)*cd(1,3)-g(15)*cd(2,3)-g(16)*cd(3,3)))
     4 +0.5*gy*(-g(18)*cd(1,3)-g(19)*cd(2,3)-g(20)*cd(3,3))
     1 -2.*lambda*(g(1)**2+g(2)**2+g(3)**2+g(4)**2-vev**2)*g(4)
     5 -0.5*((-gw*g(23)+gy*g(24))*g(3)+gw*g(21)*g(1)-gw*g(22)*g(2))
c
c W-fluxes:
c
      r(5)=0.
c
      r(6)=laplacian(6)+gw*( 
     2   -(dfdx(10)*f(14,i,j,k)-dfdx(14)*f(10,i,j,k))
     2   -(dfdy(10)*f(15,i,j,k)-dfdy(14)*f(11,i,j,k))
     2   -(dfdz(10)*f(16,i,j,k)-dfdz(14)*f(12,i,j,k))
     3   -(f(10,i,j,k)*fs(3,1,1)-f(14,i,j,k)*fs(2,1,1))
     3   -(f(11,i,j,k)*fs(3,1,2)-f(15,i,j,k)*fs(2,1,2))
     3   -(f(12,i,j,k)*fs(3,1,3)-f(16,i,j,k)*fs(2,1,3)))
c current from Higgs:
     4   +gw*(g(1)*cd(1,4)-g(2)*cd(1,3)+g(3)*cd(1,2)-g(4)*cd(1,1))
c gauge fixing terms:
     5   -dfdx(21)-gw*(g(10)*g(23)-g(14)*g(22))
c new terms in Num. Rel. method (though we have chosen
c g(5)=0=g(9)=g(13) so these terms don't play a role here).:
     6   -gw*(g(9)*fs(3,0,1)-g(13)*fs(2,0,1))    
c
      r(7)=laplacian(7)+gw*( 
     2   -(dfdx(11)*f(14,i,j,k)-dfdx(15)*f(10,i,j,k))
     2   -(dfdy(11)*f(15,i,j,k)-dfdy(15)*f(11,i,j,k))
     2   -(dfdz(11)*f(16,i,j,k)-dfdz(15)*f(12,i,j,k))
     3   -(f(10,i,j,k)*fs(3,2,1)-f(14,i,j,k)*fs(2,2,1))
     3   -(f(11,i,j,k)*fs(3,2,2)-f(15,i,j,k)*fs(2,2,2))
     3   -(f(12,i,j,k)*fs(3,2,3)-f(16,i,j,k)*fs(2,2,3)))
     4   +gw*(g(1)*cd(2,4)-g(2)*cd(2,3)+g(3)*cd(2,2)-g(4)*cd(2,1))
     5   -dfdy(21)-gw*(g(11)*g(23)-g(15)*g(22))
     6   -gw*(g(9)*fs(3,0,2)-g(13)*fs(2,0,2))    
c
      r(8)=laplacian(8)+gw*( 
     2   -(dfdx(12)*f(14,i,j,k)-dfdx(16)*f(10,i,j,k))
     2   -(dfdy(12)*f(15,i,j,k)-dfdy(16)*f(11,i,j,k))
     2   -(dfdz(12)*f(16,i,j,k)-dfdz(16)*f(12,i,j,k))
     3   -(f(10,i,j,k)*fs(3,3,1)-f(14,i,j,k)*fs(2,3,1))
     3   -(f(11,i,j,k)*fs(3,3,2)-f(15,i,j,k)*fs(2,3,2))
     3   -(f(12,i,j,k)*fs(3,3,3)-f(16,i,j,k)*fs(2,3,3)))
     4   +gw*(g(1)*cd(3,4)-g(2)*cd(3,3)+g(3)*cd(3,2)-g(4)*cd(3,1))
     5   -dfdz(21)-gw*(g(12)*g(23)-g(16)*g(22))
     6   -gw*(g(9)*fs(3,0,3)-g(13)*fs(2,0,3))    
c
      r(9)=0.
c
      r(10)=laplacian(10)+gw*( 
     2   -(dfdx(14)*f(6,i,j,k)-dfdx(6)*f(14,i,j,k))
     2   -(dfdy(14)*f(7,i,j,k)-dfdy(6)*f(15,i,j,k))
     2   -(dfdz(14)*f(8,i,j,k)-dfdz(6)*f(16,i,j,k))
     3   -(f(14,i,j,k)*fs(1,1,1)-f(6,i,j,k)*fs(3,1,1))
     3   -(f(15,i,j,k)*fs(1,1,2)-f(7,i,j,k)*fs(3,1,2))
     3   -(f(16,i,j,k)*fs(1,1,3)-f(8,i,j,k)*fs(3,1,3)))
     4   +gw*(-g(1)*cd(1,3)-g(2)*cd(1,4)+g(3)*cd(1,1)+g(4)*cd(1,2))
     5   -dfdx(22)-gw*(g(14)*g(21)-g(6)*g(23))
     6   -gw*(g(13)*fs(1,0,1)-g(5)*fs(3,0,1))    
c
      r(11)=laplacian(11)+gw*( 
     2   -(dfdx(15)*f(6,i,j,k)-dfdx(7)*f(14,i,j,k))
     2   -(dfdy(15)*f(7,i,j,k)-dfdy(7)*f(15,i,j,k))
     2   -(dfdz(15)*f(8,i,j,k)-dfdz(7)*f(16,i,j,k))
     3   -(f(14,i,j,k)*fs(1,2,1)-f(6,i,j,k)*fs(3,2,1))
     3   -(f(15,i,j,k)*fs(1,2,2)-f(7,i,j,k)*fs(3,2,2))
     3   -(f(16,i,j,k)*fs(1,2,3)-f(8,i,j,k)*fs(3,2,3)))
     4   +gw*(-g(1)*cd(2,3)-g(2)*cd(2,4)+g(3)*cd(2,1)+g(4)*cd(2,2))
     5   -dfdy(22)-gw*(g(15)*g(21)-g(7)*g(23))
     6   -gw*(g(13)*fs(1,0,2)-g(5)*fs(3,0,2))    
c
      r(12)=laplacian(12)+gw*( 
     2   -(dfdx(16)*f(6,i,j,k)-dfdx(8)*f(14,i,j,k))
     2   -(dfdy(16)*f(7,i,j,k)-dfdy(8)*f(15,i,j,k))
     2   -(dfdz(16)*f(8,i,j,k)-dfdz(8)*f(16,i,j,k))
     3   -(f(14,i,j,k)*fs(1,3,1)-f(6,i,j,k)*fs(3,3,1))
     3   -(f(15,i,j,k)*fs(1,3,2)-f(7,i,j,k)*fs(3,3,2))
     3   -(f(16,i,j,k)*fs(1,3,3)-f(8,i,j,k)*fs(3,3,3)))
     4   +gw*(-g(1)*cd(3,3)-g(2)*cd(3,4)+g(3)*cd(3,1)+g(4)*cd(3,2))
     5   -dfdz(22)-gw*(g(16)*g(21)-g(8)*g(23))
     6   -gw*(g(13)*fs(1,0,3)-g(5)*fs(3,0,3))    
c
      r(13)=0.
c
      r(14)=laplacian(14)+gw*( 
     2   -(dfdx(6)*f(10,i,j,k)-dfdx(10)*f(6,i,j,k))
     2   -(dfdy(6)*f(11,i,j,k)-dfdy(10)*f(7,i,j,k))
     2   -(dfdz(6)*f(12,i,j,k)-dfdz(10)*f(8,i,j,k))
     3   -(f(6,i,j,k)*fs(2,1,1)-f(10,i,j,k)*fs(1,1,1))
     3   -(f(7,i,j,k)*fs(2,1,2)-f(11,i,j,k)*fs(1,1,2))
     3   -(f(8,i,j,k)*fs(2,1,3)-f(12,i,j,k)*fs(1,1,3)))
     4   +gw*(g(1)*cd(1,2)-g(2)*cd(1,1)-g(3)*cd(1,4)+g(4)*cd(1,3))
     5   -dfdx(23)-gw*(g(6)*g(22)-g(10)*g(21))
     6   -gw*(g(5)*fs(2,0,1)-g(9)*fs(1,0,1))    
c
      r(15)=laplacian(15)+gw*( 
     2   -(dfdx(7)*f(10,i,j,k)-dfdx(11)*f(6,i,j,k))
     2   -(dfdy(7)*f(11,i,j,k)-dfdy(11)*f(7,i,j,k))
     2   -(dfdz(7)*f(12,i,j,k)-dfdz(11)*f(8,i,j,k))
     3   -(f(6,i,j,k)*fs(2,2,1)-f(10,i,j,k)*fs(1,2,1))
     3   -(f(7,i,j,k)*fs(2,2,2)-f(11,i,j,k)*fs(1,2,2))
     3   -(f(8,i,j,k)*fs(2,2,3)-f(12,i,j,k)*fs(1,2,3)))
     4   +gw*(g(1)*cd(2,2)-g(2)*cd(2,1)-g(3)*cd(2,4)+g(4)*cd(2,3))
     5   -dfdy(23)-gw*(g(7)*g(22)-g(11)*g(21))
     6   -gw*(g(5)*fs(2,0,2)-g(9)*fs(1,0,2))    
c
      r(16)=laplacian(16)+gw*( 
     2   -(dfdx(8)*f(10,i,j,k)-dfdx(12)*f(6,i,j,k))
     2   -(dfdy(8)*f(11,i,j,k)-dfdy(12)*f(7,i,j,k))
     2   -(dfdz(8)*f(12,i,j,k)-dfdz(12)*f(8,i,j,k))
     3   -(f(6,i,j,k)*fs(2,3,1)-f(10,i,j,k)*fs(1,3,1))
     3   -(f(7,i,j,k)*fs(2,3,2)-f(11,i,j,k)*fs(1,3,2))
     3   -(f(8,i,j,k)*fs(2,3,3)-f(12,i,j,k)*fs(1,3,3)))
     4   +gw*(g(1)*cd(3,2)-g(2)*cd(3,1)-g(3)*cd(3,4)+g(4)*cd(3,3))
     5   -dfdz(23)-gw*(g(8)*g(22)-g(12)*g(21))
     6   -gw*(g(5)*fs(2,0,3)-g(9)*fs(1,0,3))    
c
c Y-fluxes:
c
        r(17)=0.
c
        r(18)=laplacian(18)
     1   +gy*(g(1)*cd(1,2)-g(2)*cd(1,1)+g(3)*cd(1,4)-g(4)*cd(1,3))
     5   -dfdx(24)
c
        r(19)=laplacian(19)
     1   +gy*(g(1)*cd(2,2)-g(2)*cd(2,1)+g(3)*cd(2,4)-g(4)*cd(2,3))
     5   -dfdy(24)
c
        r(20)=laplacian(20)
     1   +gy*(g(1)*cd(3,2)-g(2)*cd(3,1)+g(3)*cd(3,4)-g(4)*cd(3,3))
     5   -dfdz(24)
c
! Gauge condition fixing fluxes are set to 0
      r(21)=0.
      r(22)=0.
      r(23)=0.
      r(24)=0.

! Flux for |\Phi|
        tem_mod_flx=lap_mod_phi+0.5*(
     1             -conjg(d2xhatn(1,i,j,k))*cmplx(g(1),g(2))
     1             -conjg(d2xhatn(2,i,j,k))*cmplx(g(3),g(4))
     1             -conjg(d2yhatn(1,i,j,k))*cmplx(g(1),g(2))
     1             -conjg(d2yhatn(2,i,j,k))*cmplx(g(3),g(4))
     1             -conjg(d2zhatn(1,i,j,k))*cmplx(g(1),g(2))
     1             -conjg(d2zhatn(2,i,j,k))*cmplx(g(3),g(4))
     2             -conjg(cmplx(g(1),g(2)))*d2xhatn(1,i,j,k)
     2             -conjg(cmplx(g(3),g(4)))*d2xhatn(2,i,j,k)
     2             -conjg(cmplx(g(1),g(2)))*d2yhatn(1,i,j,k)
     2             -conjg(cmplx(g(3),g(4)))*d2yhatn(2,i,j,k)
     2             -conjg(cmplx(g(1),g(2)))*d2zhatn(1,i,j,k)
     2             -conjg(cmplx(g(3),g(4)))*d2zhatn(2,i,j,k)
     3             -2.*(conjg(dxhatn(1,i,j,k)))*cmplx(dfdx(1),dfdx(2))
     3             -2.*(conjg(dxhatn(2,i,j,k)))*cmplx(dfdx(3),dfdx(4))
     3             -2.*(conjg(dyhatn(1,i,j,k)))*cmplx(dfdy(1),dfdy(2))
     3             -2.*(conjg(dyhatn(2,i,j,k)))*cmplx(dfdy(3),dfdy(4))
     3             -2.*(conjg(dzhatn(1,i,j,k)))*cmplx(dfdz(1),dfdz(2))
     3             -2.*(conjg(dzhatn(2,i,j,k)))*cmplx(dfdz(3),dfdz(4))
     4             -2.*conjg(cmplx(dfdx(1),dfdx(2)))*dxhatn(1,i,j,k)
     4             -2.*conjg(cmplx(dfdx(3),dfdx(4)))*dxhatn(2,i,j,k)
     4             -2.*conjg(cmplx(dfdy(1),dfdy(2)))*dyhatn(1,i,j,k)
     4             -2.*conjg(cmplx(dfdy(3),dfdy(4)))*dyhatn(2,i,j,k)
     4             -2.*conjg(cmplx(dfdz(1),dfdz(2)))*dzhatn(1,i,j,k)
     4             -2.*conjg(cmplx(dfdz(3),dfdz(4)))*dzhatn(2,i,j,k)
     5             +conjg(hatn(1,i,j,k))*cmplx(
     5              r(1)-laplacian(1),r(2)-laplacian(2))
     5             +conjg(hatn(2,i,j,k))*cmplx(
     5              r(3)-laplacian(3),r(4)-laplacian(4))
     6             +conjg(cmplx(r(1)-laplacian(1),r(2)-laplacian(2)))
     6             *hatn(1,i,j,k)
     6             +conjg(cmplx(r(3)-laplacian(3),r(4)-laplacian(4)))
     6             *hatn(2,i,j,k))

!Check step to ensure |\Phi| only has real components
          if(aimag(tem_mod_flx).lt.(0.000001)) then
               r(25)=tem_mod_flx
          else
               print*, '|\Phi| flux has imaginary components'
               stop
          endif

      return
      end
