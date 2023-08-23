!===========================================
!SU(2)xU(1) field strengths
!The 0 components are trivial for the static
!relaxation case 
!===========================================
c
      subroutine fieldstrength(f,i,j,k,dfdx,dfdy,dfdz,fs)
c
      implicit none
      include 'parameters.inc'
      integer i,j,k
      real*8 f,r
      real*8 dfdx,dfdy,dfdz,fs
c
      dimension f(nf,1:localsize,1:localsize,1:localsize)
      dimension dfdx(nf),dfdy(nf),dfdz(nf)
      dimension fs(4,0:3,0:3)
c
c field strength: W^a_{\mu\nu}=fs(a,mu,nu):
c
c a=1:
      fs(1,0,0)=0.
      fs(1,0,1)=0.
      fs(1,0,2)=0.
      fs(1,0,3)=0.
      fs(1,1,0)=-fs(1,0,1)
      fs(1,1,1)=0.
      fs(1,1,2)=dfdx(7)-dfdy(6)
     1         +gw*(f(10,i,j,k)*f(15,i,j,k)-f(11,i,j,k)*f(14,i,j,k))
           fs(2,1,2)=dfdx(11)-dfdy(10)
      fs(1,1,3)=dfdx(8)-dfdz(6)
     1         +gw*(f(10,i,j,k)*f(16,i,j,k)-f(12,i,j,k)*f(14,i,j,k))
      fs(1,2,0)=-fs(1,0,2)
      fs(1,2,1)=-fs(1,1,2)
      fs(1,2,2)=0.
      fs(1,2,3)=dfdy(8)-dfdz(7)
     1         +gw*(f(11,i,j,k)*f(16,i,j,k)-f(12,i,j,k)*f(15,i,j,k))
      fs(1,3,0)=-fs(1,0,3)
      fs(1,3,1)=-fs(1,1,3)
      fs(1,3,2)=-fs(1,2,3)
      fs(1,3,3)=0.
c a=2:
      fs(2,0,0)=0.
      fs(2,0,1)=0.
      fs(2,0,2)=0.
      fs(2,0,3)=0.
      fs(2,1,0)=-fs(2,0,1)
      fs(2,1,1)=0.
      fs(2,1,2)=dfdx(11)-dfdy(10)
     1         +gw*(f(14,i,j,k)*f(7,i,j,k)-f(15,i,j,k)*f(6,i,j,k))
      fs(2,1,3)=dfdx(12)-dfdz(10)
     1         +gw*(f(14,i,j,k)*f(8,i,j,k)-f(16,i,j,k)*f(6,i,j,k))
      fs(2,2,0)=-fs(2,0,2)
      fs(2,2,1)=-fs(2,1,2)
      fs(2,2,2)=0.
      fs(2,2,3)=dfdy(12)-dfdz(11)
     1         +gw*(f(15,i,j,k)*f(8,i,j,k)-f(16,i,j,k)*f(7,i,j,k))
      fs(2,3,0)=-fs(2,0,3)
      fs(2,3,1)=-fs(2,1,3)
      fs(2,3,2)=-fs(2,2,3)
      fs(2,3,3)=0.
c a=3:
      fs(3,0,0)=0.
      fs(3,0,1)=0.
      fs(3,0,2)=0.
      fs(3,0,3)=0.
      fs(3,1,0)=-fs(3,0,1)
      fs(3,1,1)=0.
      fs(3,1,2)=dfdx(15)-dfdy(14)
     1         +gw*(f(6,i,j,k)*f(11,i,j,k)-f(7,i,j,k)*f(10,i,j,k))
      fs(3,1,3)=dfdx(16)-dfdz(14)
     1         +gw*(f(6,i,j,k)*f(12,i,j,k)-f(8,i,j,k)*f(10,i,j,k))
      fs(3,2,0)=-fs(3,0,2)
      fs(3,2,1)=-fs(3,1,2)
      fs(3,2,2)=0.
      fs(3,2,3)=dfdy(16)-dfdz(15)
     1         +gw*(f(7,i,j,k)*f(12,i,j,k)-f(8,i,j,k)*f(11,i,j,k))
      fs(3,3,0)=-fs(3,0,3)
      fs(3,3,1)=-fs(3,1,3)
      fs(3,3,2)=-fs(3,2,3)
      fs(3,3,3)=0.
c
c hypercharge field strength:
      fs(4,0,0)=0.
      fs(4,0,1)=0.
      fs(4,0,2)=0.
      fs(4,0,3)=0.
      fs(4,1,0)=-fs(4,0,1)
      fs(4,1,1)=0.
      fs(4,1,2)=dfdx(19)-dfdy(18)
      fs(4,1,3)=dfdx(20)-dfdz(18)
      fs(4,2,0)=-fs(4,0,2)
      fs(4,2,1)=-fs(4,1,2)
      fs(4,2,2)=0.
      fs(4,2,3)=dfdy(20)-dfdz(19)
      fs(4,3,0)=-fs(4,0,3)
      fs(4,3,1)=-fs(4,1,3)
      fs(4,3,2)=-fs(4,2,3)
      fs(4,3,3)=0.

      return
      end
