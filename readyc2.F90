!     =================
      SUBROUTINE READYC_UDEL4
!     =================
!     PREPARATION OF BAROTROPIC AND BAROCLINIC INTEGRATION
 
!     ADVECTION + DIFFUSION
 
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
use work_mod
use tracer_mod
use pmix_mod
use POP_HaloMod
use POP_GridHorzMod
#if ( defined TIDEMIX )
use forc_mod, only: su,sv,USTAR,BUOYTUR, BUOYSOL,wave_dis
#else
use forc_mod, only: su,sv,USTAR,BUOYTUR, BUOYSOL
#endif
use domain
use grid
use blocks
use advection
use operators
use LICOM_Error_mod
#ifdef BIHAR
use hmix_del4
#else
use hmix_del2
#endif
use msg_mod 
use gather_scatter
use distribution
use constant_mod
!use canuto_2010_mod
use canuto_mod
      IMPLICIT NONE
!      REAL(r8)  :: WKP (KMP1)
      INTEGER   :: IWK,n2, iblock,kmb,iwk1 !LPF20160729
      REAL(r8)  :: WK1 (KM-1) ,WK2 (KM-1), WK3 (KM-1),WK4(KM)
      REAL(r8)  :: WP1 (KM) ,WP2 (KM), WP3 (KM) !LPF20160728
      !REAL(r8)  :: WP1 (KM-1) ,WP2 (KM-1), WP3 (KM-1) !LPF20160728
      !REAL(r8)  :: WP1 (KM) ,WP2 (KM), WP3 (KM)
      REAL(r8)  :: WP4 (KM) ,WP5 (KM), WP6 (KM)
      !REAL(r8)  :: WP4 (KMM1) ,WP5 (KMM1), WP6 (KMM1)
      REAL(r8)  :: WP7 (KM) ,WP8(KM),tau_mag,mldtmp !LPF20160715
      !REAL(r8)  :: WP7 (KMM1) ,WP8(KM),tau_mag !LPF20160715
      REAL(r8)  :: WP9 ,WP10, WP11
      REAL(r8),dimension(IMT,JMT,KM,MAX_BLOCKS_CLINIC) :: WP12,WP13,riv1,riv2
      REAL(r8)  :: epsln,RKV,RKV1,ek0
      REAL(r8)  :: adv_x1,adv_x2,adv_x,adv_y1,adv_y2,adv_z,diff_u1,diff_u2,diff_v1,diff_v2
      REAL(r8)  :: dlux,dlvx,dluy,dlvy,dluz,dlvz,adv_z1,adv_z2,adv_z3,adv_z4
      real(r8)  :: hdvk(imt,jmt), hduk(imt,jmt), adv_uu(imt,jmt,km), adv_vv(imt,jmt,km)
      real(r8)  :: tq,sq
!LPF20160505
#ifdef BCKMEX
      real(r8) :: diff_back(imt,jmt,max_blocks_clinic),&
                  diff_back_sh(imt,jmt,max_blocks_clinic),&
                  diff_back_nh(imt,jmt,max_blocks_clinic) 
#endif
      REAL(r8)    :: DENS
      EXTERNAL DENS

!
!     real (r8) :: ttt(imt_global, jmt_global)
      type (block) :: this_block          ! block information for current block
      integer ::  ErrorCode

       
#if (defined CANUTO)
      REAL(r8)  :: AIDIF
#endif

!---------------------------------------------------------------------
!     COMPUTE THE HORIZONTAL VISCOSITY
!---------------------------------------------------------------------
 
#if (defined BIHAR)

!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block,K,HDUK,HDVK)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      this_block = get_block(blocks_clinic(iblock),iblock)
      DO K = 1,KM
          call hdiffu_del4(k, HDUK, HDVK, up(:,:,k,iblock), vp(:,:,k,iblock), this_block)
          do j = 3, jmt-2
          do i = 3, imt-2
             dlv (i,j,k,iblock) = dlv(i,j,k,iblock) + hdvk(i,j)
             dlu (i,j,k,iblock) = dlu(i,j,k,iblock) + hduk(i,j)
          end do
          end do
      END DO
   END DO

#else
!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block,hduk,hdvk,K)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      this_block = get_block(blocks_clinic(iblock),iblock)
      DO K = 1,KM
         call hdiffu_del2(k, HDUK, HDVK, up(:,:,k,iblock), vp(:,:,k,iblock), this_block)
         DO J = 3, JMT-2
            DO I = 3,IMT-2
               dlv (i,j,k,iblock) = dlv(i,j,k,iblock) + hdvk(i,j)
               dlu (i,j,k,iblock) = dlu(i,j,k,iblock) + hduk(i,j)
            END DO
         END DO
      END DO
   END DO
 
#endif
 
      RETURN
      END SUBROUTINE READYC_UDEL4
 
 
