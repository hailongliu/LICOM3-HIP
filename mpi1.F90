subroutine mpi1(E1,E2)

#include <def-undef.h>
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

real (r8) :: E1 ,E2

call  mpi_reduce(E1,E2,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)

end subroutine mpi1
