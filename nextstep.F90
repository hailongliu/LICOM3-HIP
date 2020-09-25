!     ============================================
      subroutine NEXTSTEP 
#include <def-undef.h>
!     ============================================
!     move by Lin Pengfei & Liu Hailong 2013 October
use param_mod
use pconst_mod
use dyn_mod
use tracer_mod
use forc_mod
use work_mod
use domain

implicit none
!
     integer :: iblock
!
!$OMP PARALLEL DO PRIVATE (iblock,k,j,i)
   do iblock = 1, nblocks_clinic
      DO k = 1,km
         DO j = 1,jmt ! Dec. 5, 2002, Yongqiang Yu
            DO i = 1,imt
               up (i,j,k,iblock) = u (i,j,k,iblock)
               vp (i,j,k,iblock) = v (i,j,k,iblock)
               utf (i,j,k,iblock) = u (i,j,k,iblock)
               vtf (i,j,k,iblock) = v (i,j,k,iblock)
               atb (i,j,k,1,iblock) = at (i,j,k,1,iblock)
               atb (i,j,k,2,iblock) = at (i,j,k,2,iblock)
            END DO
         END DO
      END DO
   END DO

         CALL VINTEG (U,UB)
         CALL VINTEG (V,VB)

!$OMP PARALLEL DO PRIVATE (iblock,j,i)
    do iblock = 1, nblocks_clinic
      DO j = 1,jmt ! Dec. 5, 2002, Yongqiang Yu
         DO i = 1,imt
            h0p (i,j,iblock)= h0 (i,j,iblock)
            ubp (i,j,iblock)= ub (i,j,iblock)
            vbp (i,j,iblock)= vb (i,j,iblock)
            h0f (i,j,iblock)= h0 (i,j,iblock)
            h0bf (i,j,iblock)= h0 (i,j,iblock)
         END DO
      END DO
   end do
 
     ISB = 0
     ISC = 0
     IST = 0
! 
     return
      END SUBROUTINE NEXTSTEP
