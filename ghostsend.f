        subroutine ghostsend(f,pid,np,bb,fb,lb,rb,db,ub,itime)

!=======================================================================
!              Subroutine to exchange ghostpoints
! Note: nf is explicity entered as an integer in mpi_send commands
!=======================================================================

                use mpi
                implicit none
                include 'parameters.inc'
                
        integer ghostsize,ierr,np,status(mpi_status_size)
        integer bb,fb,lb,rb,db,ub,itime
        integer pid,pidb,pidf,pidl,pidr,pidd,pidu,proci,procj,prock
        integer msgsend

        real*8 f
        dimension f(nf,1:localsize,1:localsize,1:localsize)
        

        call processcoord(pid,proci,procj,prock)

! Following are the previous and next sublattice pids for a given pid
! in x,y,z directions.

        pidb=pid-1
        pidf=pid+1

        pidl=pid-nprocessx
        pidr=pid+nprocessx

        pidd=pid-(nprocessx*nprocessx)
        pidu=pid+(nprocessx*nprocessx)

        msgsend=itime

        ghostsize=nf*3*(localsize-12)*(localsize-12)

!------------------------Send Data in X-direction---------------------------------
!Send backwards and forwards everywhere except boundaries
        if(proci/=0.and.proci/=nprocessx-1) then
                call mpi_send(f(1:25,fb-3:fb-1,lb:rb,db:ub),
     1                  ghostsize,mpi_double_precision,pidf,
     2                  msgsend,mpi_comm_world,ierr)
                call mpi_send(f(1:25,bb+1:bb+3,lb:rb,db:ub),
     1                  ghostsize,mpi_double_precision,pidb,
     2                  msgsend,mpi_comm_world,ierr)
        
!Send only forwards at the oth sublattice
        elseif(proci==0) then
                call mpi_send(f(1:25,fb-3:fb-1,lb:rb,db:ub),
     1                  ghostsize,mpi_double_precision,pidf,
     2                  msgsend,mpi_comm_world,ierr)

!Send only backwards at the last(nprocessx^3th) sublattice
        else
                call mpi_send(f(1:25,bb+1:bb+3,lb:rb,db:ub),
     1                  ghostsize,mpi_double_precision,pidb,
     2                  msgsend,mpi_comm_world,ierr)
        endif

!------------------------Send Data in Y-direction---------------------------------

        if(procj/=0.and.procj/=nprocessx-1) then
                call mpi_send(f(1:25,bb:fb,rb-3:rb-1,db:ub),
     1                  ghostsize,mpi_double_precision,pidr,
     2                  msgsend,mpi_comm_world,ierr)
                call mpi_send(f(1:25,bb:fb,lb+1:lb+3,db:ub),
     1                  ghostsize,mpi_double_precision,pidl,
     2                  msgsend,mpi_comm_world,ierr)
        elseif(procj==0) then
                call mpi_send(f(1:25,bb:fb,rb-3:rb-1,db:ub),
     1                  ghostsize,mpi_double_precision,pidr,
     2                  msgsend,mpi_comm_world,ierr)

        else
                call mpi_send(f(1:25,bb:fb,lb+1:lb+3,db:ub),
     1                  ghostsize,mpi_double_precision,pidl,
     2                  msgsend,mpi_comm_world,ierr)
        endif

!------------------------Send Data in Z-direction---------------------------------

        if(prock/=0.and.prock/=nprocessx-1) then
                call mpi_send(f(1:25,bb:fb,lb:rb,ub-3:ub-1),
     1                  ghostsize,mpi_double_precision,pidu,
     2                  msgsend,mpi_comm_world,ierr)
                call mpi_send(f(1:25,bb:fb,lb:rb,db+1:db+3),
     1                  ghostsize,mpi_double_precision,pidd,
     2                  msgsend,mpi_comm_world,ierr)
        elseif(prock==0) then
                call mpi_send(f(1:25,bb:fb,lb:rb,ub-3:ub-1),
     1                  ghostsize,mpi_double_precision,pidu,
     2                  msgsend,mpi_comm_world,ierr)

        else
                call mpi_send(f(1:25,bb:fb,lb:rb,db+1:db+3),
     1                  ghostsize,mpi_double_precision,pidd,
     2                  msgsend,mpi_comm_world,ierr)
        endif

        return
        end
