!  CVS: $Id: k2_3.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
#include <def-undef.h>
 
#if (defined ISO)
!     ===============
      SUBROUTINE K2_3
!     ===============
 
!     compute "K2(,,3)" at the center of the northern face of "T" cells
!     use "c1e10" to keep the exponents in range.
use precision_mod 
use param_mod
use isopyc_mod
use constant_mod, only : PI
use pconst_mod
use domain
use grid
      IMPLICIT NONE
 
      REAL(r8):: c1e10,eps,p5,c0,c1,p25,fxd,fxe,fxc,fxa,fxb,chkslp,olmask,xx
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
      c0 = 0.0D0
      c1 = 1.0D0
      p25 = 0.25D0
 
!yyq 0803408
!-----------------------------------------------------------------------
!     d(rho_bary_barz)/dz centered on northern face of "T" cells
!     Note: values involving ocean surface and ocean bottom are
!           estimated afterwards using a linear extrapolation
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (iblock,j,k,m,fxd)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO k = 2,km -1
            m = kisrpl (k)
            fxd = c1e10* p25* dzr (k)
            DO i = 2,imt-1
               e (i,k,j,3,iblock) = fxd * (rhoi (i,k -1,j,m,iblock) - rhoi (i,k +1,j,m,iblock) &
                             + rhoi (i,k -1,j +1,m,iblock) - rhoi (i,k +1,j +1,m,iblock))
            END DO
         END DO
      END DO
  end do
 
!-----------------------------------------------------------------------
!     linearly extrapolate densities to ocean surface for calculation
!     of d(rho_bary_barz)/dz involving level 1.
!-----------------------------------------------------------------------
 
!     k   = 1
      fxd = c1e10* dzr (1)
      fxe = dzw (0) + dzw (1)
      m = kisrpl (1)
!$OMP PARALLEL DO PRIVATE (iblock,j,fxa,fxb,fxc)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO i = 2,imt-1
            fxa = p5* (rhoi (i,2,j,m,iblock) + rhoi (i,2,j +1,m,iblock))
            fxb = p5* (rhoi (i,1,j,m,iblock) + rhoi (i,1,j +1,m,iblock))
            fxc = dzwr (1)* (fxb * fxe- fxa * dzw (0))
            e (i,1,j,3,iblock) = fxd * (fxc - p5* (fxa + fxb))
         END DO
      END DO
   end do
 
 
!-----------------------------------------------------------------------
!     linearly extrapolate densities to ocean bottom for calculation
!     of d(rho_bary_barz)/dz involving bottom level.
!-----------------------------------------------------------------------
 
      e(:,km,:,3,:) = c0
 
 
!$OMP PARALLEL DO PRIVATE (iblock,j,k,m,fxa,fxb,fxc,fxe)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO i = 2,imt-1
            k = min (kmt (i,j,iblock),kmt (i,j +1,iblock))
            IF (k /= 0) THEN
               fxe = dzw (k -1) + dzw (k)
               m = kisrpl (k)
               fxa = p5* (rhoi (i,k -1,j,m,iblock) + rhoi (i,k -1,j +1,m,iblock))
               fxb = p5* (rhoi (i,k,j,m,iblock) + rhoi (i,k,j +1,m,iblock))
               fxc = dzwr (k -1)* (fxb * fxe- fxa * dzw (k))
               e (i,k,j,3,iblock) = dzr (k)* c1e10* (p5* (fxa + fxb) - fxc)
            END IF
         END DO
      END DO
  end do 
 
!-----------------------------------------------------------------------
!     "e(,,,1)" = d(rho_barx_bary)/dx centered on north face of "T" cell
!     "e(,,,2)" = d(rho)/dy on north face of "T" cells
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (iblock,j,k,m)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO k = 1,km
            m = kisrpl (k)
            DO i = 2,imt-1
               e (i,k,j,1,iblock) = p5*c1e10/(dxu(i,j,iblock)+dxu(i+1,j,iblock))*(rhoi (i+1,k,j+1,m,iblock) -      &
                                    rhoi (i -1,k,j +1,m,iblock)+ rhoi (i +1,k,j,m,iblock) - rhoi(i-1,k,j,m,iblock))
               e (i,k,j,2,iblock) = tmask (i,k,j,iblock)* tmask (i,k,j +1,iblock)/hue(i,j,iblock)  &
                            * c1e10 * (rhoi (i,k,j +1,m,iblock) - rhoi (i,k,j,m,iblock))  
            END DO
         END DO
      END DO
   end do 
 
!-----------------------------------------------------------------------
!     if any one of the 4 neighboring corner grid points is a land point
!     set "e(i,k,j,1)" to zero. note that "e(i,k,j,1)" will be used
!     only in the slope check.
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (iblock,j,k,olmask)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO k = 1,km
            DO i = 2,imt-1
               olmask = tmask (i -1,k,j +1,iblock)* tmask (i +1,k,j +1,iblock)        &
                       * tmask (i -1,k,j,iblock) * tmask (i +1,k,j,iblock) 
               if (olmask < c1) e (i,k,j,1,iblock) = c0
            END DO
         END DO
      END DO
   end do 
 
!-----------------------------------------------------------------------
!     impose zonal boundary conditions at "i"=1 and "imt"
!-----------------------------------------------------------------------
 
!-----------------------------------------------------------------------
!     compute "K2", using "slmxr" to limit vertical slope of isopycnal
!     to guard against numerical instabilities.
!-----------------------------------------------------------------------
 
#ifdef LDD97
!lhl060506
!$OMP PARALLEL DO PRIVATE (iblock,j,k,chkslp,SLOPEMOD,NONDIMR)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO k = 1,km
            DO i = 1,imt-1
               chkslp = - sqrt (e (i,k,j,1,iblock)**2+ e (i,k,j,2,iblock)**2)* slmxr
               if (e (i,k,j,3,iblock) > chkslp) e (i,k,j,3,iblock) = chkslp
!
               SLOPEMOD= sqrt(e(i,k,j,1,iblock)**2+e(i,k,j,2,iblock)**2)/abs(e(i,k,j,3,iblock)+eps)
               F1(i,j,k,iblock)=0.5D0*( 1.0D0 + tanh((0.004D0-SLOPEMOD)/0.001D0))
               NONDIMR=-ZKT(k)/(RRD2(i,j,iblock)*(SLOPEMOD+eps))
!
               IF ( NONDIMR>=1.0 ) THEN
               F2(i,j,k,iblock)=1.0D0
               ELSE
               F2(i,j,k,iblock) = 0.5D0*( 1.0D0 + SIN(PI*(NONDIMR-0.5D0)))
               ENDIF
               K2 (i,k,j,3,iblock) = ( - e (i,k,j,2,iblock)* e (i,k,j,3,iblock)* &
                                F1(i,j,k,iblock)*F2(i,j,k,iblock) ) &
                              / (e (i,k,j,3,iblock)**2+ eps)
            END DO
         END DO
      END DO
   end do
!lhl060506

#else

!$OMP PARALLEL DO PRIVATE (j,k,chkslp)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO k = 1,km
            DO i = 1,imt-1
               chkslp = - sqrt (e (i,k,j,1,iblock)**2+ e (i,k,j,2,iblock)**2)* slmxr
               if (e (i,k,j,3,iblock) > chkslp) e (i,k,j,3,iblock) = chkslp
               K2 (i,k,j,3,iblock) = ( - e (i,k,j,2,iblock)* e (i,k,j,3,iblock)* fzisop (k)) &
                              / (e (i,k,j,3,iblock)**2+ eps) 
            END DO
         END DO
      END DO
   end do
#endif

 
      RETURN
      END SUBROUTINE K2_3
#else
      SUBROUTINE K2_3
      RETURN
      END SUBROUTINE K2_3
#endif 
