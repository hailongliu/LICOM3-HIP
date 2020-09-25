#define LOGMSG()
      module smuvh

#include <def-undef.h>
use precision_mod
use param_mod
use msg_mod
use domain
use constant_mod
use grid
use POP_GridHorzMod
use POP_HaloMod

      IMPLICIT NONE
!
      integer            :: ErrorCode   ! temporary

      public :: smuv_3d, smts, smz0, smuv_2d


      contains

!     ========================
      SUBROUTINE SMUV_2D(X,Z,fil_lat)
!     ========================
!     1-D zonal smoother


      INTEGER :: JFS1,JFS2,JFN1,JFN2,KK, NCY
      REAL(r8)    :: fil_lat
      REAL(r8)    :: X (IMT,JMT,max_blocks_clinic),XS (IMT),Z (IMT,JMT,max_blocks_clinic)
      INTEGER     :: NN(JMT), MAX_NN, iblock
!
!
      MAX_NN = 10
      do j =3, jmt-2
         if (cos(ulat(1,j,1)).le.cos(fil_lat*DEGtoRAD)) then
            NN(j) = int(cos(fil_lat*DEGtoRAD)/abs(cos(ulat(1,j,1)))*1.2D0)
         else 
            NN(j) = 0
         endif
      enddo

!
      DO NCY = 1,MAX_NN
!
!$OMP PARALLEL DO PRIVATE (iblock,k,j)
   do iblock = 1, nblocks_clinic
         do j = 3, jmt-2
            if (NN(j) .ge. NCY) then       
                  DO I = 1,IMT
                     XS (I)= X (I,J,iblock)* Z (I,J,iblock)
                  END DO
                  DO I = 3,imt-2
                     X (I,J,iblock)= (0.5D0*XS(I)+0.25D0*(XS(I-1)+XS(I+1)))*Z(I,J,iblock)
                  END DO
            endif
         end do
   end do
!
         call POP_HaloUpdate(X, POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
!
   END DO

      RETURN
      END SUBROUTINE SMUV_2D





!     ========================
      SUBROUTINE SMUV_3D (X,Z,fil_lat)
!     ========================
!     1-D zonal smoother


      INTEGER :: JFS1,JFS2,JFN1,JFN2,KK, NCY
      REAL(r8)    :: fil_lat
      REAL(r8)    :: X (IMT,JMT,KM,max_blocks_clinic),XS (IMT),Z (IMT,JMT,KM,max_blocks_clinic)
      INTEGER     :: NN(JMT), MAX_NN, iblock

!
      MAX_NN = 10
      do j =3, jmt-2
         if (cos(ulat(1,j,1)).le.cos(fil_lat*DEGtoRAD)) then
            NN(j) = int(cos(fil_lat*DEGtoRAD)/cos(ulat(1,j,1))*1.2D0)
         else 
            NN(j) = 0
         endif
      enddo

!
      DO NCY = 1,MAX_NN
!
!$OMP PARALLEL DO PRIVATE (iblock,k,j)
   do iblock = 1, nblocks_clinic
         do j = 3, jmt-2
            if (NN(j) .ge. NCY) then       
               DO K = 1,KM
                  DO I = 1,IMT
                     XS (I)= X (I,J,K,iblock)* Z (I,J,K,iblock)
                  END DO
                  DO I = 3,imt-2
                     X (I,J,K,iblock)= (0.5D0*XS(I)+0.25D0*(XS(I-1)+XS(I+1)))*Z(I,J,K,iblock)
                  END DO
               ENDDO
            endif
         END DO
     END DO
!
         call POP_HaloUpdate(X, POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
!
   END DO

      RETURN
      END SUBROUTINE SMUV_3D



!     ========================
      SUBROUTINE SMTS (X,Z,fil_lat)
!     ========================
!     1-D zonal smoother


      INTEGER :: JFS1,JFS2,JFN1,JFN2,KK, NCY
      REAL(r8)    :: fil_lat, ek0
      REAL(r8)    :: X (IMT,JMT,KM,max_blocks_clinic),XS (IMT),Z (IMT,JMT,KM,max_blocks_clinic)
      INTEGER     :: NN(JMT), MAX_NN, iblock


!     fil_lat=66.D0
!lhl      fil_lat=56.D0
!

      MAX_NN = 10
      nn =0
      do j =3, jmt-2
         if (cos(tlat(1,j,1)).le.cos(fil_lat*DEGtoRAD)) then
            NN(j) = min(max_nn,abs(int(cos(fil_lat*DEGtoRAD)/cos(tlat(1,j,1))*1.2D0)))
         endif
         if ( nn(j) < 0) nn(j) = 0
      enddo
!!$OMP PARALLEL DO PRIVATE (iblocks,K,J,xs)
   DO NCY = 1,MAX_NN
         do iblock = 1, nblocks_clinic
         do j = 3, jmt-2
            if (NN(j) .ge. NCY) then
               DO K = 1,KM
                  DO I = 1,IMT
                     XS (I)= X (I,J,K,iblock)* Z (I,J,K,iblock)
                  END DO
                  DO I = 3, imt-2
                     X(I,J,K,iblock)=(XS(I)*(1.0D0-0.25D0*Z(I-1,J,K,iblock)-0.25D0*Z(I+1,J,K,iblock)) &
                               +0.25D0*(XS(I-1)+XS(I+1)))*Z(I,J,K,iblock)
                  END DO
               ENDDO
           end if
         end do
         end do
         call POP_HaloUpdate(X , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
  END DO


      RETURN
      END SUBROUTINE SMTS



!     ========================
      SUBROUTINE SMZ0 (X,Z,fil_lat)
!     ========================
!     1-D zonal smoother


      INTEGER :: JFS1,JFS2,JFN1,JFN2,NCY
      REAL(r8)    :: fil_lat
      REAL(r8)    :: X (IMT,JMT,max_blocks_clinic),XS (IMT),Z (IMT,JMT,max_blocks_clinic)
      INTEGER    :: NN(JMT), MAX_NN, iblock


!lhl      fil_lat=54.D0
!
         
      MAX_NN = 10
      do j = 3,jmt-2
         if (cos(tlat(1,j,1)).le.cos(fil_lat*DEGtoRAD)) then
            NN(j) = int(cos(fil_lat*DEGtoRAD)/abs(cos(tlat(1,j,1)))*1.2D0)
         else
            NN(j) = 0
         endif
      enddo

!!$OMP PARALLEL DO PRIVATE (iblock,J,nn,xs)
      DO NCY = 1,MAX_NN
         do iblock = 1, nblocks_clinic
         do j =3, jmt-2
            if (NN(j) .ge. NCY) then
               DO I = 1,IMT
                  XS (I)= X (I,J,IBLOCK)* Z (I,J,iblock)
               END DO
               DO I = 3, imt-2
                  X(I,J,iblock)=(XS(I)*(1.0D0-0.25D0*Z(I-1,J,iblock)-0.25D0*Z(I+1,J,iblock)) &
                            +0.25D0*(XS(I-1)+XS(I+1)))*Z(I,J,iblock)
               END DO
            endif
         enddo
         enddo
!        
         call POP_HaloUpdate(X , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
!
      END DO


   RETURN
   END SUBROUTINE SMZ0

!

    end module smuvh
