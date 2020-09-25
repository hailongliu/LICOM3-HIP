!  CVS: $Id: isoadv.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
#include <def-undef.h>
 
#if (defined ISO)
!     =================
      SUBROUTINE ISOADV
!     =================
 
!     compute isopycnal transport velocities.
use precision_mod 
use param_mod
use pconst_mod
use isopyc_mod
use msg_mod
use domain
use grid

      IMPLICIT NONE
 
      REAL(r8):: p5,c0,fxa
      INTEGER :: jstrt, iblock
 
!-----------------------------------------------------------------------
!     set local constants
!-----------------------------------------------------------------------
      p5 = 0.5D0
 
      c0 = 0.0D0
 
!-----------------------------------------------------------------------
!     compute the meridional component of the isopycnal mixing velocity
!     at the center of the northern face of the "t" cells.
!-----------------------------------------------------------------------
      adv_vntiso=c0 
      adv_vbtiso=c0 
      adv_vetiso=c0 
!$OMP PARALLEL DO PRIVATE (iblock,fxa)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO k = 2,km -1
            DO i = 2,imt-1
               fxa = - p5* dzr(k)*athkdf(i,j,iblock)
               adv_vntiso (i,k,j,iblock) = fxa * tmask (i,k,j,iblock)* tmask (i,k,j +1,iblock)* ( &
                                    K2 (i,k -1,j,3,iblock) - K2 (i,k +1,j,3,iblock))
            END DO
         END DO
      END DO
   END DO
 
 
!     consider the top and bottom levels. "K2" is assumed to be zero
!     at the ocean top and bottom.
 
!     k = 1
!$OMP PARALLEL DO PRIVATE (iblock,fxa)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO i = 2,imt-1
            fxa = - p5* dzr (1)* athkdf(i,j,iblock)
            adv_vntiso (i,1,j,iblock) = - fxa * tmask (i,1,j,iblock)* tmask (i,1,j +1,iblock)&
                                 * (K2 (i,1,j,3,iblock) + K2 (i,2,j,3,iblock))   
         END DO
      END DO
   END DO
 
 
 
 
!$OMP PARALLEL DO PRIVATE (iblock)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO i = 2, imt-1
            k = min (KMT(i,j,iblock),KMT(i,j +1,iblock))
            IF (k /= 0) THEN
               adv_vntiso (i,k,j,iblock) = - p5*dzr(k)*athkdf(i,j,iblock)*tmask(i,k,j,iblock) &
                                    * tmask (i,k,j +1,iblock)* (K2 (i,k,j,3,iblock)   &
                                      + K2 (i,k -1,j,3,iblock))          
            END IF
         END DO
      END DO
   END DO

 
!-----------------------------------------------------------------------
!     compute the zonal component of the isopycnal mixing velocity
!     at the center of the eastern face of the "t" grid box.
!-----------------------------------------------------------------------
 
      jstrt = 2
 
!$OMP PARALLEL DO PRIVATE (iblock,fxa)
   do iblock = 1, nblocks_clinic
      DO j = 2, jmt-1
         DO k = 2,km -1
            DO i = 2,imt-1
               fxa = - p5* dzr (k)* athkdf(i,j,iblock)
               adv_vetiso (i,k,j,iblock) = fxa * tmask (i,k,j,iblock)* tmask (i +1,k,j,iblock) &
                                    * (K1 (i,k -1,j,3,iblock) - K1 (i,k +1,j,3,iblock))      
            END DO
         END DO
      END DO
   END DO
 
 
!     consider the top and bottom levels. "K1" is assumed to be zero
!     at the ocean top and bottom.
 
!     k = 1
!$OMP PARALLEL DO PRIVATE (iblock,fxa)
   do iblock = 1, nblocks_clinic
      DO j = 2, jmt-1
         DO i = 2,imt-1
            fxa = - p5* dzr (1)* athkdf(i,j,iblock)
            adv_vetiso (i,1,j,iblock) = - fxa * tmask (i,1,j,iblock)* tmask (i +1,1,j,iblock) &
                                 * (K1 (i,1,j,3,iblock) + K1 (i,2,j,3,iblock)) 
         END DO
      END DO
   END DO
 
 
!$OMP PARALLEL DO PRIVATE (iblock)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO i = 2, imt-1
            k = min (kmt (i,j,iblock),kmt(i +1,j,iblock))
            IF (k /= 0) THEN
               adv_vetiso (i,k,j,iblock) = - p5*dzr(k)*athkdf(i,j,iblock)*tmask(i,k,j,iblock) &
                                    * tmask(i+1,k,j,iblock)* (K1 (i,k,j,3,iblock)   &
                                      + K1 (i,k -1,j,3,iblock)) 
            END IF
         END DO
      END DO
  END DO
!----------------------------------------------------------------------
!     compute the vertical component of the isopycnal mixing velocity
!     at the center of the bottom face of the "t" cells, using the
!     continuity equation for the isopycnal mixing velocities
!-----------------------------------------------------------------------
 
 
 
!$OMP PARALLEL DO PRIVATE (iblock)
   do iblock = 1, nblocks_clinic
      DO j = 3, jmt-2
         DO k = 1,km -1
            DO i = 3,imt-2
               adv_vbtiso (i,k,j,iblock) = DZP (k)* ( &
               (adv_vetiso (i,k,j,iblock)*htw(i+1,j,iblock) - adv_vetiso (i-1,k,j,iblock)*htw(i,j,iblock))* tarea_r(i,j,iblock) + &
               (adv_vntiso (i,k,j,iblock)*hts(i,j,iblock) - adv_vntiso (i,k,j -1,iblock)*hts(i,j-1,iblock))*tarea_r(i,j,iblock))     
            END DO
         END DO
      END DO
   END DO
 
 
!$OMP PARALLEL DO PRIVATE (iblock,j)
   do iblock = 1, nblocks_clinic
      DO j = 3, jmt-2
         DO k = 1,km -1
            DO i = 3,imt-2
               adv_vbtiso (i,k,j,iblock) = adv_vbtiso (i,k,j,iblock) + adv_vbtiso (i,k -1,j,iblock)
            END DO
         END DO
      END DO
  END DO
 
 
!$OMP PARALLEL DO PRIVATE (iblock,j,i)
   do iblock = 1, nblocks_clinic
      DO j = 3, jmt-2
         DO i = 3, imt-2
            adv_vbtiso (i,kmt (i,j,iblock),j,iblock) = c0
         END DO
      END DO
   END DO

      RETURN
      END SUBROUTINE ISOADV
 
 
#else
      SUBROUTINE ISOADV
      RETURN
      END SUBROUTINE ISOADV
#endif 
