      SUBROUTINE BAROTR_F
!     =================
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
use work_mod
use msg_mod
use domain
use grid
use blocks
use hmix_del2
use hmix_del4
use operators
use smuvh
use POP_GridHorzMod
use POP_HaloMod
use global_reductions
use gather_scatter
use distribution
use constant_mod
      IMPLICIT NONE

      INTEGER :: IEB,NC,IEB_LOOP
      real(r8)    :: gstar ,am_viv,fil_lat1,fil_lat2, ek0, maxz0, minz0
      integer :: iblock,ii1,jj1,ii2,jj2, irec
      real(r8):: hduk(imt,jmt) , hdvk(imt,jmt), gradx(imt,jmt),grady(imt,jmt), div_out(imt,jmt)
      real(r8):: hdtk(imt,jmt)
      type(block):: this_block

!---------------------------------------------------------------------
!     COMPUTE THE "ARTIFICIAL" HORIZONTAL VISCOSITY
!---------------------------------------------------------------------
#if (defined BIHAR)
!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block, hduk,hdvk)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         this_block = get_block(blocks_clinic(iblock),iblock)
         call hdiffu_del4(1, HDUK, HDVK, ubp(:,:,iblock), vbp(:,:,iblock), this_block)
         DO J = 3, jmt-2
            DO I = 3, imt-2
               WKA (I,J,5,IBLOCK)= hduk(i,j)
               WKA (I,J,6,IBLOCK)= hdvk(i,j)
            END DO
         END DO
     END DO
!
!!!!!!
#else
!$OMP PARALLEL DO PRIVATE (IBLOCK, this_block, hduk, hdvk)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         this_block = get_block(blocks_clinic(iblock),iblock)
         call hdiffu_del2(1, HDUK, HDVK, ubp(:,:,iblock), vbp(:,:,iblock), this_block)
         DO J = 3, jmt-2
            DO I = 3, imt-2
               WKA (I,J,5,IBLOCK)= hduk(i,j)
               WKA (I,J,6,IBLOCK)= hdvk(i,j)
            END DO
         END DO
    END DO

#endif
      RETURN
      END SUBROUTINE BAROTR_F

      SUBROUTINE BAROTR_TDEL
!     =================
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
use work_mod
use msg_mod
use domain
use grid
use blocks
use hmix_del2
use hmix_del4
use operators
use smuvh
use POP_GridHorzMod
use POP_HaloMod
use global_reductions
use gather_scatter
use distribution
use constant_mod
      IMPLICIT NONE

      INTEGER :: IEB,NC,IEB_LOOP
      real(r8)    :: gstar ,am_viv,fil_lat1,fil_lat2, ek0, maxz0, minz0
      integer :: iblock,ii1,jj1,ii2,jj2, irec
      real(r8):: hduk(imt,jmt) , hdvk(imt,jmt), gradx(imt,jmt),grady(imt,jmt), div_out(imt,jmt)
      real(r8):: hdtk(imt,jmt)
      REAL (r8)   :: DT2K(imt,jmt)
      type(block):: this_block

         if ( mod(nc,2) == 0 ) then
!2020 YYQ
#ifdef  BIHAR
!2020 YYQ
             !call hdifft_del4(1,hdtk,h0p(:,:,iblock),this_block)
             call hdifft_del4(1,dt2k,hdtk,h0p(:,:,iblock),this_block)
!2020 YYQ
#else
             call hdifft_del2(1,hdtk,h0p(:,:,iblock),this_block)
#endif
         else
             hdtk = c0
         end if
         work(:,:,1)=hdtk(:,:)
      RETURN
      END SUBROUTINE BAROTR_TDEL


