!  CVS: $Id: smag.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!      this program is Smagrinsky horizontal viscosity
!      based on MOM2 Manual and Rosati A. & K. Miyakoda (1988)
!      writted by liu hailong jun 2001.
 
#include <def-undef.h>
#if ( defined SMAG)
      SUBROUTINE smag2 (KK)
use precision_mod 
use param_mod
use pconst_mod
use work_mod
 
      IMPLICIT none
      INTEGER :: KK
      REAL(r8)    :: abcd,bond
 
#if (defined SMAG_FZ)
!$OMP PARALLEL DO PRIVATE (iblock,J,I)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmm
         DO i = 2,imm
            wka (i,j,13,iblock)= ((0.5D0*oux(j)*(-wka(i-1,j,11,iblock)+wka(i+1,j,11,iblock)) &
                          - (r1e (j)* wka (i,j,12,iblock) - r1f (j)* wka (i,j-1,12,iblock) &
                          +r1e(j+1)*wka(i,j+1,12,iblock)-r1f(j+1)*wka(i,j,12,iblock))))
            wka (i,j,14,iblock)= -((0.5D0*oux(j)*(-wka(i-1,j,12,iblock)+wka(i+1,j,12,iblock)) &
                          +(r1e(j)*wka(i,j,11,iblock)-r1f(j)* wka(i,j-1,11,iblock)  &
                          + r1e(j+1)*wka (i,j+1,11,iblock)-r1f(j+1)*wka(i,j,11,iblock))))
         END DO
      END DO
   END DO
!
!-new
#else
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmm
         DO i = 2,imm
            wka (i,j,13,iblock)= (0.5D0* oux (j)* (wka (i +1,j,11,iblock) - wka (i -1,j,11,iblock)) &
                        - r1e(j)*wka(i,j +1,12,iblock) + r1f (j)* wka (i,j -1,12,iblock))
            wka (i,j,14,iblock)= - (0.5D0* oux (j)*(wka(i+1,j,12,iblock)-wka(i-1,j,12,iblock)) &
                          + r1e (j)* wka (i,j +1,11,iblock) - r1f (j)* wka (i,j -1,11,iblock))
         END DO
      END DO
   END DO
!
#endif
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   do iblock = 1, nblocks_clinic
      DO j = jst,jmt
         DO i = 1,imt
            wka (i,j,15,iblock)= sqrt (2D0*wka (i,j,13,iblock)**2+ 2D0*wka (i,j,14,iblock)**2)
         END DO
      END DO
   END DO
 
 
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   do iblock = 1, nblocks_clinic
      DO j = 1,jmt
         DO i = 1,imt
            amx (i,j,KK,iblock)=  wka (i,j,15,iblock)* cxu (j)       
            amy (i,j,KK,iblock)=  wka (i,j,15,iblock)* cyu (j)       
         END DO
      END DO
   END DO

 
 
!$OMP PARALLEL DO PRIVATE (J,I)
   do iblock = 1, nblocks_clinic
      DO j = jst,jmt
         DO i = 1,imt
            wka (i,j, 7,iblock)= amx (i,j,KK,iblock)* wka (i,j,13,iblock)* viv (i,j,KK,iblock)
            wka (i,j, 8,iblock)= amy (i,j,KK,iblock)* wka (i,j,14,iblock)* viv (i,j,KK,iblock)
            wka (i,j, 9,iblock)= amx (i,j,KK,iblock)* wka (i,j,14,iblock)* viv (i,j,KK,iblock)
            wka (i,j,10,iblock)= - amy (i,j,KK,iblock)* wka (i,j,13,iblock)* viv (i,j,KK,iblock)
         END DO
      END DO
   END DO
 
      END SUBROUTINE smag2
 
 
 
      SUBROUTINE smag3
use precision_mod 
use param_mod
use pconst_mod
use work_mod
use dyn_mod
      IMPLICIT none
 
      REAL(r8)    :: dst (imt,jmt,2,iblock),dd (imt,jmt,iblock),abcd
 
!$OMP PARALLEL DO PRIVATE (iblock)
     do iblock = 1, nblocks_clinic
      DO k = 1,km
#if (defined SMAG_FZ)
         DO j = 2,jmt
            DO i = 2,imm
               dst(i,j,1,iblock)=0.5D0*otx(j)*((utl(i,j,k,iblock)-utl(i-1,j,k,iblock))+ & 
                                 (utl(i+1,j,k,iblock)-utl(i,j,k,iblock))) &
                         -(r1e(j)*vtl(i,j,k,iblock)-r1f(j)*vtl(i,j-1,k,iblock)+r1e(j)*vtl(i,j+1,k,iblock)       &
                         -r1f(j+1)*vtl(i,j,k,iblock))
               dst(i,j,2,iblock)=-(0.5D0*otx(j)*((vtl(i,j,k,iblock)-vtl(i-1,j,k,iblock))+ & 
                            (vtl(i+1,j,k,iblock)-vtl(i,j,k,iblock))) &
                          +(r1e(j)*utl(i,j,k,iblock)-r1f(j)*utl(i,j-1,k,iblock)+r1e(j+1)*utl(i,j+1,k,iblock)      &
                          -r1f(j+1)*utl(i,j,k,iblock)))
            END DO
         END DO
#else
!
         DO j = 2,jmm
            DO i = 2,imm
               dst (i,j,1,iblock)= viv(i,j,k,iblock)*(0.5D0*oux(j)*(utl(i+1,j,k,iblock)-utl(i-1,j,k,iblock)) &
                           -r1e(j)*vtl(i,j+1,k,iblock)+r1f(j)*vtl(i,j-1,k,iblock))     
               dst (i,j,2,iblock)= - viv(i,j,k,iblock)*(0.5D0*oux(j)*(vtl(i+1,j,k,iblock)-vtl(i-1,j,k,iblock)) &
                            +r1e (j)* utl (i,j +1,k,iblock) - r1f (j)* utl (i,j -1,k,iblock))   
            END DO
         END DO
!
 
#endif
 
         DO j = jsm,jmm
            DO i = 1,imt
               dd (i,j,iblock)= sqrt (dst (i,j,1,iblock)**2+ dst (i,j,2,iblock)**2)
            END DO
         END DO
 
 
         DO j = jst,jmt
            DO i = 1,imt
               am3 (i,j,k,iblock)= am (j)
            END DO
         END DO
 
 
         DO j = jsm,jmm
            DO i = 1,imt
               am3 (i,j,k,iblock)= am3 (i,j,k,iblock)*rr+dd (i,j,iblock)* cxu(j)*rr    
            END DO
         END DO
!
      END DO
 
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
   do iblock = 1, nblocks_clinic
      DO k = 1,km
         DO j = jsm,jmm
            DO i = 2,imm
               IF (vit (i,j,k,iblock) > 0.5) THEN
                  ah3 (i,j,k,iblock)= vit (i,j,k,iblock)* (am3 (i,j,k,iblock)* viv (i,j,k,iblock) &
                             + am3 (i +1,j,k,iblock)* viv (i +1,j,k,iblock) &
                             + am3 (i,j -1,k,iblock)* viv (i,j -1,k,iblock) &
                             + am3 (i +1,j -1,k,iblock)* viv (i +1,j -1,k,iblock))/ &
                             (viv (i,j,k,iblock) + viv (i +1,j,k,iblock) &
                             + viv (i,j -1,k,iblock) + viv (i +1,j -1,k,iblock))          
               END IF
            END DO
         END DO
      END DO
  END DO
!
 
      END SUBROUTINE smag3
 
 
#else
      SUBROUTINE smag2 ()
      END SUBROUTINE smag2
 
 
      SUBROUTINE smag3 ()
      END SUBROUTINE smag3
 
#endif 
