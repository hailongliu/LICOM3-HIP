
subroutine mpi2()

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use tracer_mod
use dyn_mod
use work_mod
use pmix_mod
use msg_mod
use forc_mod, only: psa,USTAR,BUOYTUR, BUOYSOL,NSWV,SWV
use domain
use grid
use blocks
use constant_mod
use operators
implicit none

call mpi_barrier(mpi_comm_ocn,ierr)
end subroutine mpi2
