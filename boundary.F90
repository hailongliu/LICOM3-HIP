!  CVS: $Id: boundary.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      subroutine boundary
!     ================
!     To compute bounary of subdomain for the each processor.
 
#include <def-undef.h>
use param_mod
use pconst_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn, &
                   nbrs,my_rank,UPUP,DOWN,RIGHT,LEFT,dims,coords,cartcomm,reorder,periods
use dyn_mod, only : buffer
use domain
use grid
use POP_HaloMod
use POP_GridHorzMod
use grid
      IMPLICIT NONE
    integer :: iblock, errorCode
!
      vit = 0.0_r8
      viv = 0.0_r8
!
   call POP_HaloUpdate(KMT, POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode,fillValue = 0_i4)
!
!$OMP PARALLEL DO PRIVATE (IBLOCK)
   do iblock = 1,nblocks_clinic
   do j=1,jmt-1
   do i=2,imt
      KMU(i,j,iblock) = min(KMT(i-1,j  ,iblock),KMT(i,j  ,iblock), &
                            KMT(i-1,j+1,iblock),KMT(i,j+1,iblock))
   end do
   end do
   end do

   call POP_HaloUpdate(KMU, POP_haloClinic, POP_gridHorzLocSWCorner,&
                       POP_fieldKindScalar, errorCode,fillValue = 0_i4)
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,j,i)
      do iblock = 1, nblocks_clinic
      do j=  1, jmt
      do i=  1, imt
         do k = 1, kmt(i,j,iblock)
            vit(i,j,k,iblock) = 1.0_r8
         end do
         do k = 1, kmu(i,j,iblock)
            viv(i,j,k,iblock) = 1.0_r8
         end do
      end do
      end do
      end do
!

      end subroutine boundary
 
 
