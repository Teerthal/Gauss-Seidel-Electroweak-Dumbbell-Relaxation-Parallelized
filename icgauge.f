!================Initial conditions for the gauge constraints===============!

        subroutine icgamma(f,bx,fx,ly,ry,dz,uz,lxb,lyl,lzd)

        implicit none
        include 'parameters.inc'
        
        integer i,j,k,bx,fx,ly,ry,dz,uz,lxb,lyl,lzd

        real*8 f
        real*8 dfdx,dfdy,dfdz,d2fdx2,d2fdy2,d2fdz2

        dimension f(nf,1:localsize,1:localsize,1:localsize)
        dimension dfdx(nf),dfdy(nf),dfdz(nf)
        dimension d2fdx2(nf),d2fdy2(nf),d2fdz2(nf)

        do k=dz,uz
            do j=ly,ry
                do i=bx,fx

        call derivatives(f,i,j,k,dfdx,dfdy,dfdz,
     1                  d2fdx2,d2fdy2,d2fdz2,lxb,lyl,lzd)

                f(21,i,j,k)=dfdx(6)+dfdy(7)+dfdz(8)
                f(22,i,j,k)=dfdx(10)+dfdy(11)+dfdz(12)
                f(23,i,j,k)=dfdx(14)+dfdy(15)+dfdz(16)
                f(24,i,j,k)=dfdx(18)+dfdy(19)+dfdz(20)
                enddo
            enddo
        enddo


        return
        end
