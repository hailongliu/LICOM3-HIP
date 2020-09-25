!  CVS: $Id: grids.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ================
      SUBROUTINE GRIDS
!     ================
!     TOPOGRAPHY & GRIDS
!-----------------------------------------------------------------------
!
! Purpose: Set up some constants dependent on model grids.
!
! Author: Yongqiang Yu and Hailong Liu, Dec. 31, 2002
!
!
!-----------------------------------------------------------------------


#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use pmix_mod
use work_mod
use constant_mod
use cdf_mod, only : start1,count1,start2,count2,start3,count3,start4,count4
#ifdef SPMD
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif
use grid
      IMPLICIT NONE
#include <netcdf.inc>

!
      REAL(r8)    :: rpart,efold1,efold2,swarg1,swarg2
      REAL(r8)    :: AJQ,rscl1,rscl2
      INTEGER :: iblock
#if ( defined TIDEMIX )
       REAL(r8)    :: fun_up(imt,jmt,max_blocks_clinic)
#endif


      if (mytid==0)then
      write(6,*)"Beginning------GRIDS !"
      endif

#if (defined SOLAR)
      rpart = 0.58D0
      efold1 = 0.35D0
      efold2 = 23.0D0
      rscl1 = 1.0D0/ efold1
      rscl2 = 1.0D0/ efold2

      DO k = 1,kmm1
         swarg1 = max (ZKP (k +1)* rscl1, -70.0)
         swarg2 = max (ZKP (k +1)* rscl2, -70.0)
         pen (k) = rpart * exp (swarg1) + (1.0D0- rpart)* exp (swarg2)
      END DO


      DO k = 1,kmm1
         pen (k) = pen (k)* OD0CP
      END DO

#endif

#if ( defined TIDEMIX )
       DO iblock = 1,max_blocks_clinic
       DO J = 1,JMT
         DO I = 1,IMT
            FUN_UP(I,J,iblock)=exp(-(abs(ZKP(int(kmt(i,j,iblock))+1))-abs(ZKP(1)))/decay_scale)/ODZT(1)
!            FUN_UP(I,J,iblock)=0.0     !yuzp-2016/11/29
            if ((vit(i,j,1,iblock).gt.0.5)) then
            DO K = 1,int(kmt(i,j,iblock))-1
            fun_up(i,j,iblock) = fun_up(i,j,iblock)+&
                    exp(-(abs(ZKP(int(kmt(i,j,iblock))+1))-abs(ZKP(K+1)))/decay_scale)/ODZT(K+1)
!                    exp(-(abs(ZKP(int(kmt(i,j,iblock))+1))-1./ODZT(1)-abs(ZKP(K+1)))/decay_scale)/ODZT(K+1)
!            fun_up(i,j,iblock) = fun_up(i,j,iblock)+exp(-(-ZKP(int(kmt(i,j,iblock))+1)+ZKP(K))/decay_scale)/ODZT(K)
!            fun_up(i,j,iblock) = fun_up(i,j,iblock)+exp(-(-ZKP(int(kmt(i,j,iblock))+1)+ZKT(K))/decay_scale)*DZP(K)
            END DO
            endif

            DO K = 1,KM
            FZ_TIDE(I,J,K,iblock)=0.0
            END DO
            if ((vit(i,j,1,iblock).gt.0.5)) then
            DO K = 1,int(kmt(i,j,iblock))-1
            FZ_TIDE(I,J,K,iblock)=exp(-(abs(ZKP(int(kmt(i,j,iblock))+1))-abs(ZKP(K+1)))/decay_scale)/fun_up(i,j,iblock)     !yuzp-2016/11/29
!            FZ_TIDE(I,J,K,iblock)=exp(-(abs(ZKP(int(kmt(i,j,iblock))+1))-abs(ZKP(K+1)))/decay_scale)/fun_up(i,j,iblock)/ODZT(K+1)     !yuzp-2016/11/29
!            FZ_TIDE(I,J,K,iblock)=exp(-(abs(ZKP(int(kmt(i,j,iblock))+1))-abs(ZKP(K+1)))/decay_scale)/(decay_scale*(1.0-exp(-(abs(ZKP(int(kmt(i,j,iblock))+1))/decay_scale))))     !yuzp-2016/11/23

!            fun_up_out(i,j,K,iblock)=fun_up(i,j,iblock)     !yuzp-2016/11/19

!            FZ_TIDE(I,J,K,iblock) = exp(-(abs(ZKP(int(kmt(i,j,iblock))+1))-1./ODZT(1)-abs(ZKP(K+1)))/decay_scale)/fun_up(i,j,iblock)
!            FZ_TIDE(I,J,K,iblock) = exp(-(-ZKP(int(kmt(i,j,iblock))+1)-1./ODZT(1)+ZKP(K+1))/decay_scale)/fun_up(i,j,iblock)/ODZT(K+1)
!            FZ_TIDE(I,J,K,iblock) = exp(-(-ZKP(int(kmt(i,j,iblock))+1)+ZKT(K))/decay_scale)/fun_up(i,j,iblock)*DZP(K)
            END DO
            endif

         END DO
      END DO
      END DO

!      if (mytid.eq.master_task) print*,zkp,zkt,dzp
!      if (mytid.eq.master_task) then
!      DO J = 1,JMT
!      DO I = 1,IMT
!      write(6,*) kmt(i,j,1)
!      write(6,*) fun_up(i,j,1)
!      write(6,*) FZ_TIDE(i,j,:,1)
!      END DO

#endif

!--------------------------------------------------------------
!     compute boundary of area where vertical mixing coefficients

      if (mytid==0)then
      write(6,*)"END------------GRIDS !"
      endif
      RETURN
      END SUBROUTINE GRIDS


