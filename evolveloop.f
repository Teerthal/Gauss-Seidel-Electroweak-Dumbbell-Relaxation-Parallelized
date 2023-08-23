!=======================================================================
!Loops that call the euler evolver over all spatial points inside the 
!domain.
!The hatn and dihatn arrays are the derivatives of \hat{\Phi}
!which are held fixed throughout the evolution.
!Only the Higgs magnitude |\Phi| or nf=25 or f(25,i,j,k) and the gauge fields are
!solved for. The nf=1:4 field components are determined by
!|\Phi|\hat{\Phi}
!=======================================================================
        subroutine evolveloop(f,bb,fb,lb,rb,db,ub,lxb,lyl,lzd,
     1                  hatn,dxhatn,dyhatn,dzhatn,
     2                  d2xhatn,d2yhatn,d2zhatn,itime,pid)

        implicit none
        include 'parameters.inc'
        
        integer itime, interval
        integer i,j,k,bb,fb,lb,rb,db,ub,lxb,lyl,lzd
        integer ilocal,jlocal,klocal
        integer n, pid

        real*8 f,coeff,r
        complex*8 hatn
        complex*8 dxhatn, dyhatn, dzhatn
        complex*8 d2xhatn, d2yhatn, d2zhatn

        dimension f(nf,1:localsize,1:localsize,1:localsize)
        dimension hatn(2,1:localsize,1:localsize,1:localsize)
        dimension dxhatn(2,1:localsize,1:localsize,1:localsize)
        dimension dyhatn(2,1:localsize,1:localsize,1:localsize)
        dimension dzhatn(2,1:localsize,1:localsize,1:localsize)
        dimension d2xhatn(2,1:localsize,1:localsize,1:localsize)
        dimension d2yhatn(2,1:localsize,1:localsize,1:localsize)
        dimension d2zhatn(2,1:localsize,1:localsize,1:localsize)
        dimension r(nf)


        coeff=49./6.

! Looping over local coordinates
        do k=db,ub
            do j=lb,rb
                do i=bb,fb
! Global coordinates (suffix local)
            ilocal = i+lxb-1
            jlocal = j+lyl-1
            klocal = k+lzd-1

        call fdflux(f,i,j,k,r,lxb,lyl,lzd,
     1                  hatn,dxhatn,dyhatn,dzhatn,
     2                  d2xhatn,d2yhatn,d2zhatn)

! Omitting end boundaries
        if(abs(ilocal).ne.latx.and.abs(jlocal)
     1                          .ne.laty.and.abs(klocal).ne.latz) then

! Scalar field evolutions
! Relax the magnitude of \Phi and keep the directions fixed

        f(25,i,j,k)= f(25,i,j,k)+relaxparam*dx**2*r(25)/coeff
        f(1,i,j,k)=f(25,i,j,k)*real(hatn(1,i,j,k))
        f(2,i,j,k)=f(25,i,j,k)*aimag(hatn(1,i,j,k))
        f(3,i,j,k)=f(25,i,j,k)*real(hatn(2,i,j,k))
        f(4,i,j,k)=f(25,i,j,k)*aimag(hatn(2,i,j,k))

! Gauge field evolutions
        do n=5,nf-1
        f(n,i,j,k)= f(n,i,j,k)+relaxparam*dx**2*r(n)/coeff
        enddo
        
        endif

                enddo
            enddo
        enddo

        return
        end
