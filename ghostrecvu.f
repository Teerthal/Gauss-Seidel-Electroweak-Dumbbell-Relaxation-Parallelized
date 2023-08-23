        subroutine ghostrecvu(f,pid,np,bb,fb,lb,rb,db,ub,itime)

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
       
!-----------------------Receive Data in X-direction--------------------

        if(proci/=0) then
                call mpi_recv(f(1:25,bb-3:bb-1,lb:rb,db:ub),
     1                  ghostsize,mpi_double_precision,pidb,
     2                  msgrecv,mpi_comm_world,status,ierr)
        endif

!------------------------Receive Data in Y-direction-------------------

        if(procj/=0) then
                call mpi_recv(f(1:25,bb:fb,lb-3:lb-1,db:ub),
     1                  ghostsize,mpi_double_precision,pidl,
     2                  msgrecv,mpi_comm_world,status,ierr)
        endif

!------------------------Receive Data in Z-direction-------------------

        if(prock/=0) then
                call mpi_recv(f(1:25,bb:fb,lb:rb,db-3:db-1),
     1                  ghostsize,mpi_double_precision,pidd,
     2                  msgrecv,mpi_comm_world,status,ierr)
        endif

        return
        end
