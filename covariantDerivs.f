! Taken from su2u1.covariantderivatives by TV and verified independently
! The Subroutine computes covariant derivatives for the Higgs field \Phi in 
! su(2)u(1) theory.
      subroutine covderiv(f,i,j,k,dfdx,dfdy,dfdz,cd)
!
      implicit none
      include 'parameters.inc'
      integer i,j,k,n
      real*8 f,fd
      real*8 dfdx,dfdy,dfdz,cd
!
      dimension f(nf,1:localsize,1:localsize,1:localsize)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension cd(0:3,4)
c
c Notation: (D_\mu \Phi)^a = cd(mu,a) where
c a=1 means real part of upper component of doublet, a=2 is imaginary
c part of upper component, a=3 is real part of lower component,
c and a=4 is imaginary part of lower component..
c By convention,
c (D_\mu \Phi)^\alpha = (\partial_\mu \phi)^\alpha
c - (gw/2){\sigma^a}_{\alpha \beta}W^b_\mu \phi^\beta 
c - \gyB_\mu\phi^\alpha,
c where, a = 1,2,3 and \alpha,\beta = 1,2(the matrice indices of pauli 
c matrices )

c fields:
c     f(1)=Re[phi(1)], f(2)=Im[phi(1)],
c     f(3)=Re[phi(2)], f(4)=Im[phi(2)].
c     f(5)=W^1_0, f(6)=W^1_1, f(7)=W^1_2, f(8)=W^1_3,
c     f(9)=W^2_0, f(10)=W^2_1, f(11)=W^2_2, f(12)=W^2_3,
c     f(13)=W^3_0, f(14)=W^3_1, f(15)=W^3_2, f(16)=W^3_3,
c     f(17)=Y_0, f(18)=Y_1, f(19)=Y_2, f(20)=Y_3;
c     f(21)=\partial_\mu W^{1\mu}, f(22)=\partial_\mu W^{2\mu},
c     f(23)=\partial_\mu W^{3\mu}, f(24)=\partial_\mu Y^\mu;
c     fd(i)=time derivative of f(i)

!!Time (D_0) components: Not needed in the relaxation algorithm since it
!!is a static problem
       cd(0,1)=0.5*gw*(f(5,i,j,k)*f(4,i,j,k)
     1        -f(9,i,j,k)*f(3,i,j,k)+f(13,i,j,k)*f(2,i,j,k))
     2        +0.5*gy*f(17,i,j,k)*f(2,i,j,k)
c
       cd(0,2)=-0.5*gw*(f(5,i,j,k)*f(3,i,j,k)
     1        +f(9,i,j,k)*f(4,i,j,k)+f(13,i,j,k)*f(1,i,j,k))
     2        -0.5*gy*f(17,i,j,k)*f(1,i,j,k)
c
       cd(0,3)=0.5*gw*(f(5,i,j,k)*f(2,i,j,k)
     1        +f(9,i,j,k)*f(1,i,j,k)-f(13,i,j,k)*f(4,i,j,k))
     2        +0.5*gy*f(17,i,j,k)*f(4,i,j,k)
c
       cd(0,4)=-0.5*gw*(f(5,i,j,k)*f(1,i,j,k)
     1        -f(9,i,j,k)*f(2,i,j,k)-f(13,i,j,k)*f(3,i,j,k))
     2        -0.5*gy*f(17,i,j,k)*f(3,i,j,k)
c
!x (D_x) components
       cd(1,1)=dfdx(1)+0.5*gw*(f(6,i,j,k)*f(4,i,j,k)
     1        -f(10,i,j,k)*f(3,i,j,k)+f(14,i,j,k)*f(2,i,j,k))
     2        +0.5*gy*f(18,i,j,k)*f(2,i,j,k)
c
       cd(1,2)=dfdx(2)-0.5*gw*(f(6,i,j,k)*f(3,i,j,k)
     1        +f(10,i,j,k)*f(4,i,j,k)+f(14,i,j,k)*f(1,i,j,k))
     2        -0.5*gy*f(18,i,j,k)*f(1,i,j,k)
c
       cd(1,3)=dfdx(3)+0.5*gw*(f(6,i,j,k)*f(2,i,j,k)
     1        +f(10,i,j,k)*f(1,i,j,k)-f(14,i,j,k)*f(4,i,j,k))
     2        +0.5*gy*f(18,i,j,k)*f(4,i,j,k)
c
       cd(1,4)=dfdx(4)-0.5*gw*(f(6,i,j,k)*f(1,i,j,k)
     1        -f(10,i,j,k)*f(2,i,j,k)-f(14,i,j,k)*f(3,i,j,k))
     2        -0.5*gy*f(18,i,j,k)*f(3,i,j,k)
c
!y (D_y) components
       cd(2,1)=dfdy(1)+0.5*gw*(f(7,i,j,k)*f(4,i,j,k)
     1        -f(11,i,j,k)*f(3,i,j,k)+f(15,i,j,k)*f(2,i,j,k))
     2        +0.5*gy*f(19,i,j,k)*f(2,i,j,k)
c
       cd(2,2)=dfdy(2)-0.5*gw*(f(7,i,j,k)*f(3,i,j,k)
     1        +f(11,i,j,k)*f(4,i,j,k)+f(15,i,j,k)*f(1,i,j,k))
     2        -0.5*gy*f(19,i,j,k)*f(1,i,j,k)
c
       cd(2,3)=dfdy(3)+0.5*gw*(f(7,i,j,k)*f(2,i,j,k)
     1        +f(11,i,j,k)*f(1,i,j,k)-f(15,i,j,k)*f(4,i,j,k))
     2        +0.5*gy*f(19,i,j,k)*f(4,i,j,k)
c
       cd(2,4)=dfdy(4)-0.5*gw*(f(7,i,j,k)*f(1,i,j,k)
     1        -f(11,i,j,k)*f(2,i,j,k)-f(15,i,j,k)*f(3,i,j,k))
     2        -0.5*gy*f(19,i,j,k)*f(3,i,j,k)
c
!z (D_z) components
       cd(3,1)=dfdz(1)+0.5*gw*(f(8,i,j,k)*f(4,i,j,k)
     1        -f(12,i,j,k)*f(3,i,j,k)+f(16,i,j,k)*f(2,i,j,k))
     2        +0.5*gy*f(20,i,j,k)*f(2,i,j,k)
c
       cd(3,2)=dfdz(2)-0.5*gw*(f(8,i,j,k)*f(3,i,j,k)
     1        +f(12,i,j,k)*f(4,i,j,k)+f(16,i,j,k)*f(1,i,j,k))
     2        -0.5*gy*f(20,i,j,k)*f(1,i,j,k)
c
       cd(3,3)=dfdz(3)+0.5*gw*(f(8,i,j,k)*f(2,i,j,k)
     1        +f(12,i,j,k)*f(1,i,j,k)-f(16,i,j,k)*f(4,i,j,k))
     2        +0.5*gy*f(20,i,j,k)*f(4,i,j,k)
c
       cd(3,4)=dfdz(4)-0.5*gw*(f(8,i,j,k)*f(1,i,j,k)
     1        -f(12,i,j,k)*f(2,i,j,k)-f(16,i,j,k)*f(3,i,j,k))
     2        -0.5*gy*f(20,i,j,k)*f(3,i,j,k)
c
      return
      end
