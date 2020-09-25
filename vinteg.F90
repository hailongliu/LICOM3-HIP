!  CVS: $Id: vinteg.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ==========================
      SUBROUTINE VINTEG (WK3,WK2)
!     ==========================
use precision_mod 
use param_mod
use pconst_mod
use domain
 
      IMPLICIT NONE
      REAL(r8):: WK3 (IMT,JMT,KM,max_blocks_clinic),WK2 (IMT,JMT,max_blocks_clinic)
      integer :: iblock
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = 1, JMT
         DO I = 1,IMT
            WK2 (I,J,IBLOCK)= 0.0D0
         END DO
      END DO
   END DO
 
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,J)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
        DO J = 1,JMT
         DO I = 1,IMT
               WK2 (I,J,IBLOCK)= WK2 (I,J,IBLOCK) + DZP (K)* OHBU(I,J,IBLOCK)*  &
                                 WK3(I,J,K,IBLOCK) *VIV (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
   END DO
 
 
      RETURN
      END SUBROUTINE VINTEG
 
 
