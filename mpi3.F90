subroutine mpi3(err_norm1)

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

integer :: I1, I2
real(r8) :: err_norm1
real(r8) :: out1
!call MPI_BCAST(err_norm1, I1, MPI_PR, I2, MPI_COMM_OCN, IERR)
I1=1
I2=0
call MPI_ALLREDUCE(err_norm1, out1, I1, MPI_PR, MPI_SUM, MPI_COMM_OCN, IERR)
err_norm1=out1
end subroutine mpi3

subroutine readycF

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
use advection

implicit none

integer :: iblock
    iblock=1
   call advection_momentum(u(:,:,:,iblock),v(:,:,:,iblock),wka(:,:,:,iblock), &
                            dlu(:,:,:,iblock),dlv(:,:,:,iblock),iblock)

end subroutine readycF
