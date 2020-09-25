!     =================
      SUBROUTINE COREIDEAL_DAILY
!     =================

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use constant_mod
use forc_mod
use dyn_mod, only: u,v,buffer
use tracer_mod, only: at
use work_mod, only: work1_g
use msg_mod
use gather_scatter
use POP_GridHorzMod
use grid
use domain

      IMPLICIT NONE
#include <netcdf.inc>

    integer :: TNUM, iblock
!
    fresh = 0.0
    seaice = 0.0
    sss  = 0.0
    sst  = 0.0
    swv = 0.0
    nswv = 0.0
    su = 0.0
    sv = 0.0

   !if(mytid==0) then
   ! DO IBLOCK = 1, NBLOCKS_CLINIC
   ! do j = 1, jmt_global
   ! do i = 1, imt_global
   !    swv(i,j,iblock) = 200.0_r8*cos(tlat(i,j,iblock)/DEGtoRAD)
   !    nswv(i,j,iblock) = -150.0_r8*cos(tlat(i,j,iblock)/DEGtoRAD)
   !    su(i,j,iblock) = -0.03_r8*cos((abs(tlat(i,j,iblock)/DEGtoRAD-18.0_r8)*4.5_r8))
   ! 
   !
  ! 
  !  enddo
  !  enddo
  !  enddo   

    DO IBLOCK = 1, NBLOCKS_CLINIC
    do j = 1, jmt
    do i = 1, imt
       !sst(i,j,iblock) = 28.0_r8*cos(tlat(i,j,iblock)/DEGtoRAD)
       swv(i,j,iblock) = 200.0_r8*cos(tlat(i,j,iblock))
       nswv(i,j,iblock) = -150.0_r8*cos(tlat(i,j,iblock))
       su(i,j,iblock) = -0.1_r8*cos((abs(tlat(i,j,iblock)/DEGtoRAD)-15.0_r8)*DEGtoRAD)
       !su(i,j,iblock) = -0.03_r8*cos((abs(tlat(i,j,iblock)-pi/18.0_r8)*4.5_r8))
       !sss(i,j,iblock) = (1.5*cos((abs(tlat(i,j,iblock)-pi/6.0_r8)*1.5_r8))-0.3_r8)*0.001_r8
    end do
    end do
    end do
!
   allocate( work1_g(imt_global,jmt_global))
   call scatter_global(su,work1_g, master_task, distrb_clinic, &
                       field_loc_center, field_type_scalar)
!   call scatter_global(sst,work1_g, master_task, distrb_clinic, &
!                       field_loc_center, field_type_scalar)
!   call scatter_global(sss,work1_g, master_task, distrb_clinic, &
!                       field_loc_center, field_type_scalar)
!   call scatter_global(swv,work1_g, master_task, distrb_clinic, &
!                       field_loc_center, field_type_scalar)
!   call scatter_global(nswv,work1_g, master_task, distrb_clinic, &
!                       field_loc_center, field_type_scalar)
!   deallocate( work1_g)

    return
end subroutine COREIDEAL_DAILY
