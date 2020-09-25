!  CVS: $Id: isoflux.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
#include <def-undef.h>
 
#if (defined ISO)
!     ==============================
      SUBROUTINE ISOFLUX (MTRACE)
!     ==============================
 
!     isopycnal diffusive tracer fluxes are computed.
use precision_mod 
use param_mod
use pconst_mod
use tracer_mod
use isopyc_mod
use work_mod
use msg_mod
use domain
use grid
use constant_mod
 
      IMPLICIT NONE
 
      INTEGER :: mtrace, iblock
      REAL(r8):: fxa,fxb,fxc,fxe
 
      allocate (work_1(imt,jmt,km,max_blocks_clinic),work_2(imt,jmt,km,max_blocks_clinic))
      allocate (temp(imt,jmt,km,max_blocks_clinic))
      allocate (work_3(imt,jmt,0:km,max_blocks_clinic))
!-----------------------------------------------------------------------
!     set local constants
!-----------------------------------------------------------------------
      m = mtrace
 
!-----------------------------------------------------------------------
!     first compute the vertical tracer flux "temp" at the northern
!     face of "t" cells.
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (iblock,k)
   do iblock = 1, nblocks_clinic
      DO k = 2,km -1
         DO j = 2, jmt-1
            DO i = 2,imt-1
               temp (i,j,k,iblock)= p25* dzr (k)* (atb (i,j +1,k -1,m,iblock) - atb ( &
                  i,j +1,k +1,m,iblock) &
                   + atb (i,j,k -1,m,iblock) - atb (i,j, k +1,m,iblock))    
            END DO
         END DO
      END DO
   end do
 
 
!-----------------------------------------------------------------------
!     now consider the top level, assuming that the surface tracer
!     values are the same as the ones at "k"=1
!-----------------------------------------------------------------------
 
      k = 1
!$OMP PARALLEL DO PRIVATE (iblock,j,i)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO i = 2, imt-1
            temp (i,j,k,iblock) = 0.25* dzr (k)* (atb (i,j +1,k,m,iblock) - atb (i,   &
                          j +1,k +1,m,iblock)+atb (i,j,k,m,iblock) - atb (i,j, k +1,m,iblock))
         END DO
      END DO
   end do
 
 
!-----------------------------------------------------------------------
!     finally, consider the bottom level. the extrapolative estimator
!     is used to compute the tracer values at the ocean bottom.
!-----------------------------------------------------------------------
 
   temp (:,:,km,:)= c0
 
 
!$OMP PARALLEL DO PRIVATE (iblock,fxa,fxb,fxc,fxe)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO i = 2,imt-1
            k = min (kmt (i,j,iblock),kmt(i,j +1,iblock))
            IF (k /= 0) THEN
               fxe = dzw (k -1) + dzw (k)
               fxa = 0.5D0* (atb (i,j +1,k -1,m,iblock) + atb (i,j,k -1,m,iblock))
               fxb = 0.5D0* (atb (i,j +1,k,m,iblock) + atb (i,j,k,m,iblock))
               fxc = dzwr (k -1)* (fxb * fxe- fxa * dzw (k))
               temp (i,j,k,iblock) = dzr (k)* (0.5D0* (fxa + fxb) - fxc)
            END IF
         END DO
      END DO
   end do
 
 
!-----------------------------------------------------------------------
!     compute of meridional tracer flux
!     first calculate the effects of purely horizontal diffusion along
!     isopycnal, the background horizontal diffusion has been computed
!     before called this subroutine. (jxz)
!     add in the effects of the along isopycnal diffusion computed
!     using "K2" component of the tensor and apply land/sea masks
!-----------------------------------------------------------------------
 
            work_1 (:,1,:,:)=  c0
 
!$OMP PARALLEL DO PRIVATE (iblock,k,j)
   do iblock = 1, nblocks_clinic
      DO k = 1,km
         DO j = 2,jmt-1
            DO i = 2,imt-1
               work_1(i,j,k,iblock)=hts(i,j,iblock)*(ahisop(i,j,iblock)*(atb(i,j+1,k,m,iblock)- &
                                    atb(i,j,k,m,iblock))/hue(i,j,iblock)+ &
                                    ahisop(i,j,iblock)*K2(i,k,j,3,iblock)*temp(i,j,k,iblock))* &
                                    vit(i,j,k,iblock)*vit(i,j+1,k,iblock) 
!
               ddy_iso(i,j,k,m,iblock)=work_1(i,j,k,iblock)
!
            END DO
         END DO
      END DO
   end do
 
 
!-----------------------------------------------------------------------
!     compute the vertical tracer flux "temp" at the eastern
!     face of "t" cells.
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (iblock,k,j)
   do iblock = 1, nblocks_clinic
      DO k = 2,km -1
         DO j = 2,jmt-1
            DO i = 2,imt-1
               temp (i,j,k,iblock)= p25*dzr(k)*(atb(i+1,j,k-1,m,iblock) - atb ( &
                             i+1,j,k+1,m,iblock)+atb(i,j,k-1,m,iblock)-atb(i,j,k+1,m,iblock))     
            END DO
         END DO
      END DO
   end do
 
 
!-----------------------------------------------------------------------
!     now consider the top level, assuming that the surface tracer
!     values are the same as the ones at "k"=1
!-----------------------------------------------------------------------
 
      k = 1
!$OMP PARALLEL DO PRIVATE (iblock,j)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO i = 1,imt-1
            temp (i,j,k,iblock)= p25* dzr (k)* (atb (i +1,j,k,m,iblock) - atb (i +1,j,&
               k +1,m,iblock) + atb (i,j,k,m,iblock) - atb (i,j,k +1,m,iblock))  
         END DO
      END DO
   end do
 
 
!-----------------------------------------------------------------------
!     finally, consider the bottom level. the extrapolative estimator
!     is used to compute the tracer values at the ocean bottom.
!-----------------------------------------------------------------------
 
            temp (:,:,km,:) = c0
 
!$OMP PARALLEL DO PRIVATE (iblock,fxa,fxb,fxc,fxe)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO i = 2,imt-1
            k = min (kmt (i,j,iblock),kmt (i +1,j,iblock))
            IF (k /= 0) THEN
               fxe = dzw (k -1) + dzw (k)
               fxa = p5* (atb (i,j,k -1,m,iblock) + atb (i +1,j,k -1,m,iblock))
               fxb = p5* (atb (i,j,k,m,iblock) + atb (i +1,j,k,m,iblock))
               fxc = dzwr (k -1)* (fxb * fxe- fxa * dzw (k))
               temp (i,j,k,iblock) = dzr (k)* (p5* (fxa + fxb) - fxc)
            END IF
         END DO
      END DO
   end do
 
 
!-----------------------------------------------------------------------
!     compute of zonal tracer flux
!     first calculate the effects of purely horizontal diffusion along
!     isopycnal, the background horizontal diffusion has been computed
!     before called this subroutine. (jxz)
!     add in the effects of the along isopycnal diffusion computed
!     using "K1" component of the tensor and apply land/sea masks
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (iblock,k,j)
   do iblock = 1, nblocks_clinic
      DO k = 1,km
         DO j = 2,jmt-1
            DO i = 2,imt-1
               work_2 (i,j,k,iblock)= htw(i+1,j,iblock)*(ahisop(i,j,iblock)*(atb(i+1,j,k,m,iblock)- & 
                                      atb(i,j,k,m,iblock))/hun(i+1,j,iblock) + &
                                      ahisop(i,j,iblock)*K1(i,k,j,3,iblock)*  &
                                      temp(i,j,k,iblock))*vit(i+1,j,k,iblock)* vit(i,j,k,iblock) 
            END DO
         END DO
      END DO
   end do
 
 
!-----------------------------------------------------------------------
!     compute the vertical tracer flux "work_3" containing the K31
!     and K32 components which are to be solved explicitly. The K33
!     component will be treated semi-implicitly
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (iblock,k)
   do iblock = 1, nblocks_clinic
      DO k = 2,km
         DO j = 2,jmt-1
            DO i = 2,imt-1
               work_3 (i,j,k -1,iblock) = ahisop(i,j,iblock)*p25*vit(i,j,k,iblock)*(K3(i,k-1,j,1,iblock)* &
               (vit(i-1,j,k ,iblock)*(atb(i,j,k,m,iblock)-atb(i-1,j,k,m,iblock))/hun(i,j,iblock) &
               + vit(i-1,j,k-1,iblock)*(atb(i,j,k-1,m,iblock)-atb(i-1,j,k-1,m,iblock))/hun(i,j,iblock)  &
               + vit(i+1,j,k ,iblock)*(atb(i+1,j,k,m,iblock)- atb(i,j,k,m,iblock))/hun(i+1,j,iblock)  &
               + vit(i+1,j,k-1,iblock)*(atb(i+1,j,k-1,m,iblock)-atb(i,j,k-1,m,iblock))/hun(i+1,j,iblock)) + &
                K3 (i,k -1,j,2,iblock)* &
               (vit (i,j-1,k ,iblock)*(atb(i,j,k,m,iblock)-atb(i,j-1,k,m,iblock))/hue(i,j-1,iblock) &
               + vit(i,j-1,k-1,iblock)*(atb(i,j,k-1,m,iblock)-atb(i,j-1,k-1,m,iblock))/hue(i,j-1,iblock) &
               + vit(i,j+1,k ,iblock)*(atb(i,j+1,k,m,iblock)-atb (i,j,k,m,iblock))/hue(i,j,iblock) &
               + vit (i,j+1,k-1,iblock)*(atb(i,j+1,k-1,m,iblock)-atb (i,j,k-1,m,iblock))/hue(i,j,iblock)) )  
            END DO
         END DO
      END DO
   end do
 
!-----------------------------------------------------------------------
!     at ocean surface the flux is set to zero to reflect the no tracer
!     flux condition. Same condition is also imposed at ocean bottom.
!-----------------------------------------------------------------------
 
            work_3 (:,:, 0,:)= c0
            work_3 (:,:,km,:)= c0
 
 
!$OMP PARALLEL DO PRIVATE (iblock,k,j)
  do iblock = 1, nblocks_clinic
      DO k = 1,km
         DO j = 3,jmt-2
            DO i = 3, imt-2
               tf (i,j,k,iblock) = tf (i,j,k,iblock) &
               + tarea_r(i,j,iblock)*(work_1 (i,j,k,iblock) - work_1 (i,j -1,k,iblock)) &
               + tarea_r(i,j,iblock)*(work_2 (i,j,k,iblock) - work_2 (i -1,j,k,iblock)) &
               + dzr (k)*(work_3 (i,j,k -1,iblock) - work_3 (i,j,k,iblock))
!
               dx_iso(i,j,k,m,iblock)= tarea_r (i,j,iblock)*(work_2 (i,j,k,iblock) - work_2 (i -1,j,k,iblock)) 
               dy_iso(i,j,k,m,iblock)= tarea_r (i,j,iblock)*(work_1 (i,j,k,iblock) - work_1 (i,j -1,k,iblock))
               dz_iso(i,j,k,m,iblock)= dzr (k) * (work_3(i,j,k -1,iblock) - work_3 (i,j,k,iblock))
!
            END DO
         END DO
      END DO
   end do
 
!-----------------------------------------------------------------------
!     compute the meridional component of the isopycnal velocity mixing
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (iblock,k,j)
   do iblock = 1, nblocks_clinic
      DO k = 1,km
         DO j = 2,jmt-1
            DO i = 2,imt-1
               work_1 (i,j,k,iblock) = adv_vntiso (i,k,j,iblock)* (atb (i,j +1,k,m,iblock)   &
                                + atb (i,j,k,m,iblock))*hts(i,j,iblock)
!
            END DO
         END DO
      END DO
   end do
 
 
!-----------------------------------------------------------------------
!     compute the meridional component of the isopycnal velocity mixing
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (iblock,k,j)
   do iblock = 1, nblocks_clinic
      DO k = 1,km
         DO j = 2,jmt-1
            DO i = 2,imt-1
               work_2 (i,j,k,iblock) = adv_vetiso (i,k,j,iblock)* (atb (i +1,j,k,m,iblock)   &
                                + atb (i,j,k,m,iblock))*htw(i+1,j,iblock)
            END DO
         END DO
      END DO
   end do
 
 
!-----------------------------------------------------------------------
!     compute the vertical component of the isopycnal velocity mixing
!-----------------------------------------------------------------------
 
            work_3 (:,:, 0,:)= c0
            work_3 (:,:,km,:)= c0
 
 
!$OMP PARALLEL DO PRIVATE (iblock,k,j)
   do iblock = 1, nblocks_clinic
      DO k = 2,km
         DO j = 2,jmt-1
            DO i = 2,imt-1
               work_3 (i,j,k -1,iblock)= adv_vbtiso (i,k -1,j,iblock)* (atb (i,j,k,m,iblock) &
                                + atb (i,j,k -1,m,iblock))
            END DO
         END DO
      END DO
   end do
 
 
!$OMP PARALLEL DO PRIVATE (iblock,k,j)
   do iblock = 1, nblocks_clinic
      DO k = 1,km
         DO j = 3, jmt-2
            DO i = 3, imt-2
               tf (i,j,k,iblock) = tf (i,j,k,iblock) &
               - p5*(work_1 (i,j,k,iblock) - work_1 (i,j -1,k,iblock) &
               +work_2(i,j,k,iblock)-work_2(i-1,j,k,iblock))*tarea_r(i,j,iblock) &
               - p5*dzr(k)*(work_3(i,j,k-1,iblock) - work_3 (i,j,k,iblock)) 
!
               ax_iso(i,j,k,m,iblock) = -P5*(work_2(i,j,k,iblock)-work_2(i-1,j,k,iblock))*tarea_r(i,j,iblock)
               ay_iso(i,j,k,m,iblock) = -P5*(work_1(i,j,k,iblock)-work_1(i,j-1,k,iblock))*tarea_r(i,j,iblock)
               az_iso(i,j,k,m,iblock) = -P5*dzr(k)*(work_3(i,j,k-1,iblock) - work_3 (i,j,k,iblock))
!
            END DO
         END DO
      END DO
   end do
!
!     if (mytid == 27 .and. m ==1) then
!        write(130,*) "ist=", ist
!        write(130,*) ((tf(i,j,7,1),i=30,39),j=3,12)
!     end if
!
      deallocate (work_1,work_2,work_3,temp)
 
      RETURN
      END SUBROUTINE ISOFLUX
 
 
#else
      SUBROUTINE ISOFLUX ()
      RETURN
      END SUBROUTINE ISOFLUX
#endif 
