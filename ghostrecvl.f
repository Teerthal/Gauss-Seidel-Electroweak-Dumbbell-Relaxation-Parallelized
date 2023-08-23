        subroutine ghostrecvl(f,pid,np,bb,fb,lb,rb,db,ub,itime)

!=======================================================================
!              Subroutine to exchange ghostpoints
!=======================================================================

                use mpi
                implicit none
                include 'parameters.inc'
                
        integer ghostsize,ierr,np,status(mpi_status_size)
        integer bb,fb,lb,rb,db,ub,itime
        integer pid,pidb,pidf,pidl,pidr,pidd,pidu,proci,procj,prock
        integer msgrecv

        real*8 f
        dimension f(nf,1:localsize,1:localsize,1:localsize)
        

        call processcoord(pid,proci,procj,prock)

        pidb=pid-1
        pidf=pid+1

        pidl=pid-nprocessx
        pidr=pid+nprocessx

        pidd=pid-(nprocessx*nprocessx)
        pidu=pid+(nprocessx*nprocessx)

        msgrecv=itime

        ghostsize=nf*3*(localsize-12)*(localsize-12)
       
!-----------------------Recieve Data in X-direction--------------------

        if(proci/=nprocessx-1) then

                call mpi_recv(f(1:25,fb+1:fb+3,lb:rb,db:ub),
     1                  ghostsize,mpi_double_precision,pidf,
     2                  msgrecv,mpi_comm_world,status,ierr)

        endif

!-----------------------Recieve Data in Y-direction---------------------

        if(procj/=nprocessx-1) then

                call mpi_recv(f(1:25,bb:fb,rb+1:rb+3,db:ub),
     1                  ghostsize,mpi_double_precision,pidr,
     2                  msgrecv,mpi_comm_world,status,ierr)

        endif

!-----------------------Recieve Data in Z-direction---------------------

        if(prock/=nprocessx-1) then

                call mpi_recv(f(1:25,bb:fb,lb:rb,ub+1:ub+3),
     1                  ghostsize,mpi_double_precision,pidu,
     2                  msgrecv,mpi_comm_world,status,ierr)

        endif

        return
        end
