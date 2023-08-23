!---This subroutine collects energy from all the processes to 
!---pid=1 via mpi_recv and adds it up

        subroutine printenergy(energy,pid,np)
                use mpi

        implicit none

        real*8 energy(13), energy0(13)
        integer i,msg,ierr,np,pid,status(mpi_status_size)
        
        character (len=200) filename4
        character (len=200) filename5
        character (len=200) filename6
        character (len=200) filename7


        include 'parameters.inc'

        msg=1

!---Every process other than pid=0 sends their sublattice energies to
!---pid=0
        if(pid/=0) then
              call mpi_send(energy(1),13,mpi_double_precision,0,msg,
     1          mpi_comm_world, ierr)
        endif

!---pid=0 receives all the energies from all non-zero processes
        if(pid==0) then
              do i=1,np-1

              call mpi_recv(energy0(1),13,mpi_double_precision,i,msg,
     1             mpi_comm_world,status, ierr)
       
              energy=energy+energy0
              enddo
        endif

!---Employ mpi barrier to wait for all the energies to sent and recieved
              call mpi_barrier(mpi_comm_world,ierr)
              call mpi_bcast(energy(1),13,mpi_double_precision,0,
     1             mpi_comm_world,ierr)
              call mpi_barrier(mpi_comm_world,ierr)


!---0th process prints all the energies
                
        if(pid==0) then
!-----energy(1) needs to be divided by np since the quantity
!-----is added over all the processes        
        print*,' =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* '
        print*,' time:           ', energy(1)/np
        print*,' =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* '
        print*,' KE:             ', energy(2)
        print*,' GE:             ', energy(3)
        print*,' EEW:            ', energy(4)
        print*,' MEW:            ', energy(5)
        print*,' EEY:            ', energy(6)
        print*,' MEY:            ', energy(7)
        print*,' PEtotalE:       ', energy(8)
        print*,' totalE:         ', energy(9)
        print*,' sphere radius:  ', energy(10)/np
        print*,' energyInSphere: ', energy(11)
        print*,' energyBoundary: ', energy(12)
        print*,' String Energy:  ', energy(13)/(2.*senlen*dx)
        print*,' '
        
!--------Write total energy at every iteration----------!

        write(filename4,"(A16)") "total_energy.dat"
        open(unit=20,file=trim(adjustl(filename4)),position='append')

        write(20,*) energy(1)/np,energy(9)
        close(20)

        endif

        return
        end
