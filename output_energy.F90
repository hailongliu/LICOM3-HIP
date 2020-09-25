subroutine output_energy(MONTH1,EK0,EA0,EB0,ET0,ES0)
use precision_mod
use param_mod
use pconst_mod
use tracer_mod
use dyn_mod
use work_mod
use grid
use constant_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
use shr_sys_mod
use domain
use blocks
use constant_mod
        implicit none
        integer::month1
        real(8)::ek0,ea0,eb0,et0,es0  

      if (mytid==0 )then
        WRITE (16,FMT='(I5,I3,6D25.15)') MONTH,IDAY,EK0,EA0,EB0,ET0,ES0
!      call flush_(6)
      end if

end subroutine output_energy
