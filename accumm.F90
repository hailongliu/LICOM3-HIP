!  CVS: $Id: accumm.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE ACCUMM
!     =================
 
#include <def-undef.h>
use param_mod
use dyn_mod
use tracer_mod
use output_mod
use pconst_mod
use isopyc_mod
#if ( defined TIDEMIX )
use forc_mod, only: su,sv,lthf,sshf,lwv,swv,fresh,runoff,seaice,wave_dis
#else
use forc_mod, only: su,sv,lthf,sshf,lwv,swv,fresh,runoff,seaice
#endif
use domain
use grid, only: FCOR, FCORT 
      IMPLICIT NONE
      integer :: iblock
      real :: Num_op1 

       if(dts_accum) then
        Num_op1=float(Num_output)
        !if(mytid==0) write(*,*) 'accum thermal step Num_op1',Num_op1
       else
        Num_op1=float(imd)
       endif

#if (defined LOWRES)
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   do iblock = 1, nblocks_clinic
      DO J = 1,JMT
         DO I = 1,IMT

           Z0MON (I,J,iblock)= Z0MON (I,J,iblock) + H0(I,J,iblock)/Num_op1
           sumon (I,J,iblock)= sumon (I,J,iblock) + su(I,J,iblock)/Num_op1        !U windstress
           svmon (I,J,iblock)= svmon (I,J,iblock) + sv(I,J,iblock)/Num_op1        !V windstress
           lthfmon (I,J,iblock)= lthfmon (I,J,iblock) + lthf (I,J,iblock)/Num_op1 !latent flux
           sshfmon (I,J,iblock)= sshfmon (I,J,iblock) + sshf (I,J,iblock)/Num_op1 !sensible flux
           lwvmon (I,J,iblock)= lwvmon (I,J,iblock) + lwv (I,J,iblock)/Num_op1    !long wave flux
           swvmon (I,J,iblock)= swvmon (I,J,iblock) + swv (I,J,iblock)/Num_op1    !shortwave flux
           freshmon(I,J,iblock)= freshmon (I,J,iblock) + fresh (I,J,iblock)/34.7*1.0e3*DZP(1)/OD0/Num_op1 !virtual slat flux
           runoffmon (I,J,iblock)= runoffmon (I,J,iblock) + runoff (I,J,iblock)/Num_op1    !shortwave flux
           mldmon (I,J,iblock)= mldmon (I,J,iblock) + amld (I,J,iblock)/Num_op1
           ifracmon (I,J,iblock)= ifracmon (I,J,iblock) + seaice(I,J,iblock)/Num_op1

           netmon (I,J,1,iblock)= netmon (I,J,1,iblock) + net(I,J,1,iblock)/Num_op1 !K/s
           netmon (I,J,2,iblock)= netmon (I,J,2,iblock) + net(I,J,2,iblock)/Num_op1 !psu/s
#if ( defined TIDEMIX )
              wavedismon (I,J,iblock)= wavedismon (I,J,iblock) + wave_dis(I,J,iblock)/Num_op1
#endif
#if ( defined CANUTOMIXOUT )
              wp10_canutomon (I,J,iblock)= wp10_canutomon (I,J,iblock) + wp10_canuto(I,J,iblock)/Num_op1     !yuzp-2016/12/8
              wp11_canutomon (I,J,iblock)= wp11_canutomon (I,J,iblock) + wp11_canuto(I,J,iblock)/Num_op1     !yuzp-2016/12/8
              fcor_canutomon (I,J,iblock)= fcor_canutomon (I,J,iblock) + fcor(I,J,iblock)/Num_op1     !yuzp-2016/12/4
              fcort_canutomon (I,J,iblock)= fcort_canutomon (I,J,iblock) + fcort(I,J,iblock)/Num_op1  !yuzp-2016/12/4
#endif

      DO K = 1,KM          !yuzp-2016/11/13
#if ( defined TIDEMIX )
              richardsonmon (I,J,K,iblock)= richardsonmon (I,J,K,iblock) + richardson(I,J,K,IBLOCK)/Num_op1     !yuzp-2016/11/13
              fztidalmon (I,J,K,iblock)= fztidalmon (I,J,K,iblock) + fztidal(I,J,K,iblock)/Num_op1     !yuzp-2016/11/13
              wp3_tidalmon (I,J,K,iblock)= wp3_tidalmon (I,J,K,iblock) + wp3_tidal(I,J,K,iblock)/Num_op1     !yuzp-2016/11/13
              ak_tide1mon (I,J,K,iblock)= ak_tide1mon (I,J,K,iblock) + ak_tide1(I,J,K,iblock)/Num_op1    !yuzp-2016/11/19
#endif

#if ( defined CANUTOMIXOUT )
              wp1_canutomon (I,J,K,iblock)= wp1_canutomon (I,J,K,iblock) + wp1_canuto(I,J,K,iblock)/Num_op1     !yuzp-2016/12/4
              wp2_canutomon (I,J,K,iblock)= wp2_canutomon (I,J,K,iblock) + wp2_canuto(I,J,K,iblock)/Num_op1     !yuzp-2016/12/4
              wp3_canutomon (I,J,K,iblock)= wp3_canutomon (I,J,K,iblock) + wp3_canuto(I,J,K,iblock)/Num_op1     !yuzp-2016/12/4
              wp4_canutomon (I,J,K,iblock)= wp4_canutomon (I,J,K,iblock) + wp4_canuto(I,J,K,iblock)/Num_op1     !yuzp-2016/12/4
              wp5_canutomon (I,J,K,iblock)= wp5_canutomon (I,J,K,iblock) + wp5_canuto(I,J,K,iblock)/Num_op1     !yuzp-2016/12/4
              wp6_canutomon (I,J,K,iblock)= wp6_canutomon (I,J,K,iblock) + wp6_canuto(I,J,K,iblock)/Num_op1     !yuzp-2016/12/4
              wp7_canutomon (I,J,K,iblock)= wp7_canutomon (I,J,K,iblock) + wp7_canuto(I,J,K,iblock)/Num_op1     !yuzp-2016/12/4
              wp8_canutomon (I,J,K,iblock)= wp8_canutomon (I,J,K,iblock) + wp8_canuto(I,J,K,iblock)/Num_op1     !yuzp-2016/12/4
              wp12_canutomon (I,J,K,iblock)= wp12_canutomon (I,J,K,iblock) + wp12_canuto(I,J,K,iblock)/Num_op1     !yuzp-2016/12/4
              wp13_canutomon (I,J,K,iblock)= wp13_canutomon (I,J,K,iblock) + wp13_canuto(I,J,K,iblock)/Num_op1     !yuzp-2016/12/4

              wk1_canutomon (I,J,K,iblock)= wk1_canutomon (I,J,K,iblock) + wk1_canuto(I,J,K,iblock)/Num_op1     !yuzp-2016/12/4
              wk2_canutomon (I,J,K,iblock)= wk2_canutomon (I,J,K,iblock) + wk2_canuto(I,J,K,iblock)/Num_op1     !yuzp-2016/12/4
              wk3_canutomon (I,J,K,iblock)= wk3_canutomon (I,J,K,iblock) + wk3_canuto(I,J,K,iblock)/Num_op1     !yuzp-2016/12/4
              wk4_canutomon (I,J,K,iblock)= wk4_canutomon (I,J,K,iblock) + wk4_canuto(I,J,K,iblock)/Num_op1     !yuzp-2016/12/4
              !fcor_canutomon (I,J,iblock)= fcor_canutomon (I,J,iblock) + fcor_canuto(I,J,iblock)/Num_op1     !yuzp-2016/12/4
              !fcort_canutomon (I,J,iblock)= fcort_canutomon (I,J,iblock) + fcort_canuto(I,J,iblock)/Num_op1  !yuzp-2016/12/4
              alpha_canutomon (I,J,K,iblock)= alpha_canutomon (I,J,K,iblock) + alpha_canuto(I,J,K,iblock)/Num_op1    !yuzp-2016/12/4
              beta_canutomon (I,J,K,iblock)= beta_canutomon (I,J,K,iblock) + beta_canuto(I,J,K,iblock)/Num_op1    !yuzp-2016/12/4
#endif

      END DO     !yuzp-2016/11/13

!
!           netmon (I,J,1,iblock)= netmon (I,J,1,iblock) + net(I,J,1,iblock)/OD0CP*DZP(1)/Num_op1 !W/m^2
!           netmon (I,J,2,iblock)= netmon (I,J,2,iblock) + net(I,J,2,iblock)/34.7*1.0e3*DZP(1)/OD0/Num_op1 !kg/m^s/s
         END DO
      END DO
   end do
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
   do iblock = 1, nblocks_clinic
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
              WSMON (I,J,K,iblock)= WSMON (I,J,K,iblock) + WS(I,J,K,iblock)/Num_op1
              TSMON (I,J,K,iblock)= TSMON (I,J,K,iblock) + AT (I,J,K,1,iblock)*vit(i,j,k,iblock)/Num_op1
              SSMON (I,J,K,iblock)= SSMON (I,J,K,iblock) + (AT (I,J,K,2,iblock)*1000.+35.)*vit(i,j,k,iblock)/Num_op1
              USMON (I,J,K,iblock)= USMON (I,J,K,iblock) + U (I,J,K,iblock)/Num_op1
              VSMON (I,J,K,iblock)= VSMON (I,J,K,iblock) - V(I,J,K,iblock)/Num_op1
              akmmon (I,J,K,iblock)= akmmon (I,J,K,iblock) +akmu(I,J,K,iblock)/Num_op1
              aktmon (I,J,K,iblock)= aktmon (I,J,K,iblock) +akt(I,J,K,1,iblock)/Num_op1
              aksmon (I,J,K,iblock)= aksmon (I,J,K,iblock) +akt(I,J,K,2,iblock)/Num_op1
#if ( defined TIDEMIX )
              aktidemon (I,J,K,iblock)= aktidemon (I,J,K,iblock) + ak_tide(I,J,K,iblock)/Num_op1
#endif
#if (defined ISO_TYPE_BF)
              athkdfmon (I,J,K,iblock)= athkdfmon (I,J,K,iblock) + athkdf(I,J,K,iblock)/Num_op1
#endif
#if (defined SMAG_OUT)
              AM3MON (I,J,K,iblock)= AM3MON (I,J,K,iblock) + AM3(I,J,K,iblock)/Num_op1
#endif
              penmon (I,J,K,iblock)= penmon (I,J,K,iblock)+penetrate(i,j,k,iblock)/Num_op1 !penetration
#ifdef ISO
              !wsmon_iso (I,J,K,IBLOCK)= wsmon_iso (I,J,K,IBLOCK)+adv_vbtiso(i,k-1,j,IBLOCK)*odzp(k)/Num_op1
              !vsmon_iso (I,J,K,IBLOCK)= vsmon_iso (I,J,K,IBLOCK)+&
              !                       0.5D0*(adv_vntiso(i,k,j,IBLOCK)+adv_vntiso(i,k,j-1,iblock))/Num_op1

#if ( defined ISOOUT )
              VNTISOMON (I,J,K,iblock)= VNTISOMON (I,J,K,iblock) + VNTISO (I,J,K,iblock)/Num_op1
              VETISOMON (I,J,K,iblock)= VETISOMON (I,J,K,iblock) + VETISO (I,J,K,iblock)/Num_op1
              VBTISOMON (I,J,K,iblock)= VBTISOMON (I,J,K,iblock) + VBTISO (I,J,K,iblock)/Num_op1
              wsmon_iso(I,J,K,IBLOCK)=VBTISOMON(I,J,K,IBLOCK) !LPF20160805 !if work then can replace the wsmon_iso in other places
              vsmon_iso(I,J,K,IBLOCK)=VNTISOMON(I,J,K,IBLOCK) !LPF20161012
#endif
#endif
            END DO
         END DO
      END DO
   end do

!$OMP PARALLEL DO PRIVATE (IBLOCK,N,K,J,I)
   do iblock = 1, nblocks_clinic
      DO N = 1,NTRA
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT

               tendmon (I,J,K,N,IBLOCK)= tendmon (I,J,K,N,IBLOCK)+tend(i,j,k,n,IBLOCK)/Num_op1 !tendency
               axmon (I,J,K,N,IBLOCK)= axmon(I,J,K,N,IBLOCK)+ax(i,j,k,n,IBLOCK)/Num_op1 !zonal adv
               aymon (I,J,K,N,IBLOCK)= aymon(I,J,K,N,IBLOCK)+ay(i,j,k,n,IBLOCK)/Num_op1 !meridional adv
               azmon (I,J,K,N,IBLOCK)= azmon(I,J,K,N,IBLOCK)+az(i,j,k,n,IBLOCK)/Num_op1 !vertical adv
               dt_convmon (I,J,K,N,IBLOCK)= dt_convmon(I,J,K,N,IBLOCK)+dt_conv(i,j,k,n,IBLOCK)/Num_op1 !vertical adv


!LPF20160819
               dxmon (I,J,K,N,IBLOCK)= dxmon (I,J,K,N,IBLOCK)+dx(i,j,k,n,IBLOCK)/Num_op1 !horizontal diffusion
               !dymon (I,J,K,N,IBLOCK)= dymon (I,J,K,N,IBLOCK)+dy(i,j,k,n,IBLOCK)
               dzmon (I,J,K,N,IBLOCK)= dzmon (I,J,K,N,IBLOCK)+dz(i,j,k,n,IBLOCK)/Num_op1 !vertical diffusion
               dt_diffmon (I,J,K,N,IBLOCK)= dt_diffmon (I,J,K,N,IBLOCK)+dt_diff(i,j,k,n,IBLOCK)/Num_op1 !vertical diffusion due to dt
!LPF20160819
!
!              ddymon (I,J,K,N,IBLOCK)= ddymon (I,J,K,N,IBLOCK)+ddy(i,j,k,n,IBLOCK)
#ifdef ISO
               axmon_iso (I,J,K,N,IBLOCK)= axmon_iso (I,J,K,N,IBLOCK)+ax_iso(i,j,k,n,IBLOCK)/Num_op1 !x-horizontal advection due to eddy
               aymon_iso (I,J,K,N,IBLOCK)= aymon_iso (I,J,K,N,IBLOCK)+ay_iso(i,j,k,n,IBLOCK)/Num_op1 !y-horizontal advection due to eddy
               azmon_iso (I,J,K,N,IBLOCK)= azmon_iso (I,J,K,N,IBLOCK)+az_iso(i,j,k,n,IBLOCK)/Num_op1 !z-horizontal advection due to eddy
!LPF20160819
               dxmon_iso (I,J,K,N,IBLOCK)= dxmon_iso (I,J,K,N,IBLOCK)+dx_iso(i,j,k,n,IBLOCK)/Num_op1 !x-horizotal diffusion
               dymon_iso (I,J,K,N,IBLOCK)= dymon_iso (I,J,K,N,IBLOCK)+dy_iso(i,j,k,n,IBLOCK)/Num_op1 !y-horizontal diffusion
               dzmon_iso (I,J,K,N,IBLOCK)= dzmon_iso (I,J,K,N,IBLOCK)+dz_iso(i,j,k,n,IBLOCK)/Num_op1 !vertical diffusion due to eddy
!LPF20160819
!
!              aaymon_iso (I,J,K,N,IBLOCK)= aaymon_iso (I,J,K,N,IBLOCK)+aay_iso(i,j,k,n,IBLOCK)/Num_op1
!              ddymon_iso (I,J,K,N,IBLOCK)= ddymon_iso (I,J,K,N,IBLOCK)+ddy_iso(i,j,k,n,IBLOCK)/Num_op1
#endif
            END DO
         END DO
      END DO
      END DO
   end do

          ERR_norm2mon=ERR_norm2mon+FW_norm2/Num_op1 !global water fluxes constraint
#endif
!
      RETURN
      END SUBROUTINE ACCUMM
 
 
