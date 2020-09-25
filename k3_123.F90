!  CVS: $Id: k3_123.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
#include <def-undef.h>
 
#if (defined ISO)
!     =================
      SUBROUTINE K3_123
!     =================
 
!     compute K2(,,,1:3) at the center of the bottom face of "T" cells
!     use "c1e10" to keep the exponents in range.
use precision_mod 
use param_mod
use constant_mod, only : PI
use pconst_mod
use isopyc_mod
use domain
use grid
 
      IMPLICIT NONE
      REAL(r8) :: c1e10,eps,p5,p25,c0,c1,chkslp,ahfctr,xx
!lhl060506
      REAL(r8) , dimension(IMT,JMT,KM,max_blocks_clinic) :: F1,F2
      REAL(r8) :: NONDIMR,SLOPEMOD
      integer :: iblock

      F1=0.0D0
      F2=0.0D0

!lhl060506
 
!-----------------------------------------------------------------------
!     set local constants
!-----------------------------------------------------------------------
 
      c1e10 = 1.0d10
 
      eps = 1.0d-25
      p5 = 0.5D0
      p25 = 0.25D0
      c0 = 0.0D0
      c1 = 1.0D0
!yyq 0803408
!M
 
!$OMP PARALLEL DO PRIVATE (iblock,j,k,m)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO k = 2,km
            m = kisrpl (k)
            DO i = 2,imt-1
               e (i,k -1,j,1,iblock) =  p25* c1e10 &
               * (tmask (i -1,k -1,j,iblock)* tmask (i,k -1,j,iblock)* (rhoi (i,k -1,j,m,iblock) &
               - rhoi (i -1,k -1,j,m,iblock))/hun(i,j,iblock) &
               + tmask (i,k -1,j,iblock)* tmask (i +1,k -1,j,iblock)* (rhoi (i +1,k -1,j,m,iblock)&
               - rhoi (i,k -1,j,m,iblock))/hun(i+1,j,iblock) &
               + tmask (i -1,k,j,iblock)* tmask (i,k,j,iblock)* (rhoi (i,k,j,m,iblock) &
               - rhoi (i -1,k,j,m,iblock))/hun(i,j,iblock) &
               + tmask (i,k,j,iblock)* tmask (i +1,k,j,iblock)* (rhoi (i +1,k,j,m,iblock) &
                                - rhoi (i,k,j,m,iblock))/hun(i+1,j,iblock))                     
               e (i,k -1,j,2,iblock) =  p25* c1e10 &
               * (tmask (i,k -1,j -1,iblock)* tmask (i,k -1,j,iblock)* (rhoi (i,k -1,j,m,iblock) &
               - rhoi (i,k -1,j -1,m,iblock))/hue(i,j-1,iblock) &
               + tmask (i,k -1,j,iblock)* tmask (i,k -1,j +1,iblock)* (rhoi (i,k -1,j +1,m,iblock)&
               - rhoi (i,k -1,j,m,iblock))/hue(i,j,iblock) &
               + tmask (i,k,j -1,iblock)* tmask (i,k,j,iblock)* (rhoi (i,k,j,m,iblock) &
               - rhoi (i,k,j -1,m,iblock))/hue(i,j-1,iblock) &
               + tmask (i,k,j,iblock)* tmask (i,k,j +1,iblock)* (rhoi (i,k,j +1,m,iblock) &
                                - rhoi (i,k,j,m,iblock))/hue(i,j,iblock))
               e (i,k -1,j,3,iblock) = dzwr (k -1)* tmask (i,k -1,j,iblock)* tmask (i,k,j,iblock)* c1e10 &
                                * (rhoi (i,k -1,j,m,iblock) - rhoi (i,k,j,m,iblock))  
            END DO
         END DO
!nickbegin
!       k = km
!nickend
            e (:,km,j,1,iblock) = c0
            e (i,km,j,2,iblock) = c0
            e (:,km,j,3,iblock) = c0
      END DO
   end do
 
 
!-----------------------------------------------------------------------
!     compute "K3", using "slmxr" to limit vertical slope of isopycnal
!     to guard against numerical instabilities.
!-----------------------------------------------------------------------
 
#ifdef LDD97
!lhl060506
!$OMP PARALLEL DO PRIVATE (iblock,j,k,chkslp,SLOPEMOD,NONDIMR,ahfctr)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO k = 1,km
            DO i = 1,imt-1
               chkslp = - sqrt (e (i,k,j,1,iblock)**2+ e (i,k,j,2,iblock)**2)* slmxr
               if (e (i,k,j,3,iblock) > chkslp) e (i,k,j,3,iblock) = chkslp
!
               SLOPEMOD= sqrt(e(i,k,j,1,iblock)**2+e(i,k,j,2,iblock)**2)/abs(e(i,k,j,3,iblock)+eps)
               F1(i,j,k,iblock)=0.5D0*( 1.0D0 + tanh((0.004D0-SLOPEMOD)/0.001D0))
               NONDIMR=-ZKP(k)/(RRD1(i,j,iblock)*(SLOPEMOD+eps))
               IF ( NONDIMR>=1.0 ) THEN
               F2(i,j,k,iblock)=1.0D0
               ELSE         
               F2(i,j,k,iblock) = 0.5D0*( 1.0D0 + SIN(PI*(NONDIMR-0.5D0)))
               ENDIF
               ahfctr = 1.0D0/ (e (i,k,j,3,iblock)**2+ eps)*F1(i,j,k,iblock)*F2(i,j,k,iblock)
               K3 (i,k,j,1,iblock) = - e (i,k,j,3,iblock)* e (i,k,j,1,iblock)* ahfctr
               K3 (i,k,j,2,iblock) = - e (i,k,j,3,iblock)* e (i,k,j,2,iblock)* ahfctr
               K3 (i,k,j,3,iblock) = (e (i,k,j,1,iblock)**2+ e (i,k,j,2,iblock)**2)* ahfctr
            END DO
         END DO   
      END DO
   end do
!lhl060506

#else

!$OMP PARALLEL DO PRIVATE (iblock,j,k,chkslp,ahfctr)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO k = 1,km
            DO i = 2,imt -1
               chkslp = - sqrt (e (i,k,j,1,iblock)**2+ e (i,k,j,2,iblock)**2)* slmxr
               if (e (i,k,j,3,iblock) > chkslp) e (i,k,j,3,iblock) = chkslp
               ahfctr = 1.0/ (e (i,k,j,3,iblock)**2+ eps)
               K3 (i,k,j,1,iblock) = - e (i,k,j,3,iblock)* e (i,k,j,1,iblock)* ahfctr
               K3 (i,k,j,2,iblock) = - e (i,k,j,3,iblock)* e (i,k,j,2,iblock)* ahfctr
               K3 (i,k,j,3,iblock) = (e (i,k,j,1,iblock)**2+ e (i,k,j,2,iblock)**2)* ahfctr
            END DO
         END DO
      END DO
  end do
#endif
 
!-----------------------------------------------------------------------
!     impose zonal boundary conditions at "i"=1 and "imt"
!-----------------------------------------------------------------------
 
      RETURN
      END SUBROUTINE K3_123
#else
      SUBROUTINE K3_123
      RETURN
      END SUBROUTINE K3_123
#endif 
