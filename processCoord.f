        subroutine processcoord(pid,proci,procj,prock)
                implicit none

        include 'parameters.inc'
                
        integer pid, proci, procj, prock

        prock=int(pid/(nprocessx*nprocessx))

        procj=int((pid-nprocessx*nprocessx*prock)/nprocessx)

        proci=pid-(nprocessx*nprocessx*prock)-(nprocessx*procj)

        return
        end
