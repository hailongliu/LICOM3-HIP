!  CVS: $Id: invtri.F90,v 1.5 2003/08/12 09:06:39 lhl Exp $
#include <def-undef.h>
 
!     =================
      SUBROUTINE INVTRIT (WK,TOPBC,DCB,AIDIF,C2DTTS)
!     =================
use param_mod
use pconst_mod
use domain
use grid, only : kmt
      IMPLICIT NONE
      real(r8)  :: WK (IMT,JMT,KM,max_blocks_clinic),TOPBC (IMT,JMT,max_blocks_clinic),DCB (IMT,JMT,KM,max_blocks_clinic)
      real(r8)  :: A8 (KM),B8 (KM),C8 (KM),D8 (KM),E8 (0:KM),F8 (0:KM)
      real(r8)  :: AIDIF,C2DTTS,G0
      INTEGER :: KZ, IBLOCK
 
!$OMP PARALLEL DO PRIVATE (iblock,a8,b8,c8,d8,e8,f8,G0,KZ)
  do iblock = 1, nblocks_clinic
      JJJ: DO J = 3, jmt-2
         III: DO I = 3, imt-2
            if (kmt (I,J,iblock) > 0 ) then
            KZ = kmt (I,J,iblock)
            DO K = 2,KZ  !lyc
               A8 (K) = DCB (I,J,K -1,iblock)* ODZT (K )* ODZP (K)* C2DTTS * AIDIF
               D8 (K) = WK (I,J,K,iblock)
            ENDDO
            DO K=2,KZ-1
               C8 (K) = DCB (I,J,K,iblock)* ODZT (K+1 )* ODZP (K)* C2DTTS * AIDIF
               B8 (K) = 1.0+ A8 (K) + C8 (K)
               E8 (K -1) = 0.0
               F8 (K -1) = 0.0
            END DO
!     B. C. AT TOP
            K = 1
            A8 (K) = ODZP (K)* C2DTTS * AIDIF
            C8 (K) = DCB (I,J,K,iblock)* ODZT (K+1)* ODZP (K)* C2DTTS * AIDIF
            B8 (K) = 1.0+ C8 (K)
            D8 (K) = WK (I,J,K,iblock)
            E8 (K -1) = 0.0
            F8 (K -1) = 0.0
!     B. C. AT BOTTOM
            B8 (KZ) = 1.0+ A8 (KZ)
            C8 (KZ) = ODZP (KZ)* C2DTTS * AIDIF
            E8 (KZ) = 0.0
            F8 (KZ) = 0.0
 
!     NOW INVERT
            DO K = KZ,1, -1
               G0 = 1.0/ (B8 (K) - C8 (K)* E8 (K))
               E8 (K -1) = A8 (K)* G0
               F8 (K -1) = (D8 (K) + C8 (K)* F8 (K))* G0
            END DO
 
!     B.C. AT SURFACE
            WK (I,J,1,iblock) = (E8 (0)* TOPBC (I,J,iblock) + F8 (0))* VIT (I,J,1,iblock)
            DO K = 2,KZ
               WK (I,J,K,iblock)= (E8 (K -1)* WK (I,J,K -1,iblock) + F8 (K -1))* VIT (I,J,K,iblock)
            END DO
         end if 
         END DO III
      END DO JJJ
  end do
!
      RETURN
      END SUBROUTINE INVTRIT
 
 
!     =================
      SUBROUTINE INVTRIU (WK,TOPBC,BOMBC,DCB,AIDIF,C2DTC)
!     =================
use param_mod
use pconst_mod
use domain
use grid, only : kmu
      IMPLICIT NONE
      real(r8) :: WK(IMT,JMT,KM,max_blocks_clinic),TOPBC(IMT,JMT,max_blocks_clinic),DCB(IMT,JMT,KM,max_blocks_clinic)
      real(r8) :: A8 (KM),B8(KM),C8(KM),D8(KM),E8(0:KM),F8(0:KM),BOMBC(IMT,JMT,max_blocks_clinic)
      real(r8):: AIDIF,C2DTC,G0
      INTEGER :: KZ,IBLOCK
 
!$OMP PARALLEL DO PRIVATE (iblock,a8,b8,c8,d8,e8,f8,G0,KZ)
  do iblock = 1, nblocks_clinic
      JJJ: DO J = 3, jmt-2
         III: DO I = 3, imt-2
            IF (kmu(I,J,iblock) > 0) then
            KZ = kmu (I,J,iblock)
            DO K = 2,KZ  !lyc
               A8 (K) = DCB (I,J,K -1,iblock)* ODZT (K )* ODZP (K)* C2DTC * AIDIF
               D8 (K) = WK (I,J,K,iblock)
            ENDDO
            DO K=2,KZ-1
               C8 (K) = DCB (I,J,K,iblock)* ODZT (K+1 )* ODZP (K)* C2DTC * AIDIF
               B8 (K) = 1.0+ A8 (K) + C8 (K)
               E8 (K -1) = 0.0
               F8 (K -1) = 0.0
            END DO
!     B. C. AT TOP
            K = 1
            A8 (K) = ODZP (K)* C2DTC * AIDIF
            C8 (K) = DCB (I,J,K,iblock)* ODZT (K+1)* ODZP (K)* C2DTC * AIDIF
            B8 (K) = 1.0+ C8 (K)
            D8 (K) = WK (I,J,K,iblock)
            E8 (K -1) = 0.0
            F8 (K -1) = 0.0
!     B. C. AT BOTTOM
            B8 (KZ) = 1.0+ A8 (KZ)
            C8 (KZ) = ODZP (KZ)* C2DTC * AIDIF
            E8 (KZ) = 0.0
            F8 (KZ) = 0.0
            D8 (KZ) = WK (I,J,KZ,iblock) - BOMBC(I,J,iblock)*ODZP(KZ)*C2DTC*AIDIF
!     NOW INVERT
            DO K = KZ,1, -1
               G0 = 1.0/ (B8 (K) - C8 (K)* E8 (K))
               E8 (K -1) = A8 (K)* G0
               F8 (K -1) = (D8 (K) + C8 (K)* F8 (K))* G0
            END DO
!     B.C. AT SURFACE
            WK (I,J,1,iblock) = (E8 (0)* TOPBC (I,J,iblock) + F8 (0))* VIV (I,J,1,iblock)
            DO K = 2,KZ
               WK (I,J,K,iblock)= (E8 (K -1)* WK (I,J,K -1,iblock) + F8 (K -1))* VIV (I,J,K,iblock)
            END DO
         end if
         END DO III
      END DO JJJ
  end do
!
      RETURN
      END SUBROUTINE INVTRIU
 
