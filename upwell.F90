!  CVS: $Id: upwell.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ===============================
      SUBROUTINE UPWELL (UWK,VWK,H0WK)
!     ===============================
 
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
use work_mod
use domain
use grid
use blocks
use operators
 
      IMPLICIT NONE
      REAL(r8):: UWK (IMT,JMT,KM,max_blocks_clinic),VWK (IMT,JMT,KM,max_blocks_clinic),H0WK (IMT,JMT,max_blocks_clinic)
      type(block) :: this_block
      integer :: iblock
 
!---------------------------------------------------------------------
!     INITIALIZE WORK ARRAYS
!---------------------------------------------------------------------
 
      allocate(uk(imt,jmt,km,max_blocks_clinic),vk(imt,jmt,km,max_blocks_clinic))

!M
   work = 0.0D0
   wka  = 0.0D0
 
!$OMP PARALLEL DO PRIVATE (IBLOCK)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      call tgrid_to_ugrid(work(:,:,iblock), h0wk(:,:,iblock),iblock)
   END DO
!
 
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 1, JMT-1
            DO I = 2,IMT
               UK (I,J,K,IBLOCK)= (1.0D0+ WORK (I,J,IBLOCK)* OHBU (I,J,IBLOCK))* UWK (I,J,K,IBLOCK)
               VK (I,J,K,IBLOCK)= (1.0D0+ WORK (I,J,IBLOCK)* OHBU (I,J,IBLOCK))* VWK (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
   END DO
 
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block,K)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         this_block = get_block(blocks_clinic(iblock),iblock)
         call div(k,wka(:,:,k,iblock),uk(:,:,k,iblock),vk(:,:,k,iblock),this_block)
      END DO
   END DO
!
 
   work = 0.0D0
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 2, JMT-1
            DO I = 2,IMT-1
               WORK (I,J,IBLOCK)= WORK(I,J,IBLOCK) - DZP (K)* WKA (I,J,K,IBLOCK)* VIT (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
   END DO
 
 
!$OMP PARALLEL DO PRIVATE (IBLOCK)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 2,KM
         DO J = 2, JMT-1
            DO I = 2,IMT-1
               WS (I,J,K,IBLOCK)= VIT (I,J,K,IBLOCK)* (WS (I,J,K -1,IBLOCK) + &
               DZP (K -1)* (WORK (I,J,IBLOCK)* OHBT (I,J,IBLOCK) + WKA (I,J,K -1,IBLOCK)))
            END DO
         END DO
      END DO
  END DO
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = 2, JMT-1
         DO I = 2, IMT-1
            WORK (I,J,IBLOCK)= 1.0D0/ (1.0D0+ H0WK (I,J,IBLOCK)* OHBT (I,J,IBLOCK))
         END DO
      END DO
   END DO
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 2,KM
         DO J = 2, JMT-1
            DO I = 2,IMT-1
               WS (I,J,K,IBLOCK)= WS (I,J,K,IBLOCK)* WORK (I,J,IBLOCK)
            END DO
         END DO
      END DO
   END DO
    

 
      deallocate(uk,vk)

      RETURN
      END SUBROUTINE UPWELL
