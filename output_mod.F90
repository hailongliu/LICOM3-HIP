!  CVS: $Id: output_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module output_mod
#include <def-undef.h>
!
use precision_mod
use param_mod
!     ------------------------------------------------------------------
!     Output Arrays
!     ------------------------------------------------------------------
#if (defined DAILYACC)
      real(r4),dimension(imt,jmt,max_blocks_clinic)::z0daily
      real(r4),dimension(imt,jmt,max_blocks_clinic)::mlddaily
      real(r4),dimension(imt,jmt,max_blocks_clinic)::lthfdaily,sshfdaily,lwvdaily,swvdaily
      real(r4),dimension(imt,jmt,max_blocks_clinic)::precdaily,evapdaily,roffdaily
      real(r4),dimension(imt,jmt,max_blocks_clinic)::sudaily,svdaily
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::tsdaily,ssdaily
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::usdaily,vsdaily,wsdaily
#if (defined DAILYBUGDET)
      real(r4),allocatable,dimension(:,:,:,:,:)::tenddaily
      real(r4),allocatable,dimension(:,:,:,:,:)::dt_diffdaily,axdaily,aydaily,azdaily
      real(r4),allocatable,dimension(:,:,:,:,:)::dx_diffdaily,dxdaily,dzdaily,dt_convdaily
      !real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dxdaily,dzdaily !,dydaily
      real(r4),allocatable,dimension(:,:,:,:)::pendaily
      !real(r4),dimension(imt,jmt,km,max_blocks_clinic)::pendaily
      real(r4),allocatable,dimension(:,:,:,:,:)::dt_restoredaily
      real(r4),allocatable,dimension(:,:,:,:)::net_daily
#endif
#endif

#if (defined LOWRES)
      real(r4),dimension(imt,jmt,max_blocks_clinic)::z0mon,himon,hdmon,qicemon
      real(r4),dimension(imt,jmt,max_blocks_clinic)::lthfmon,sshfmon,lwvmon,swvmon
      real(r4),dimension(imt,jmt,max_blocks_clinic)::sumon,svmon,runoffmon, freshmon
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wsmon,tsmon,ssmon,usmon,vsmon
      real(r4),dimension(imt,jmt,2,max_blocks_clinic):: icmon
      real(r4),dimension(imt,jmt,NTRA,max_blocks_clinic)::netmon
      real(r4),dimension(imt,jmt,max_blocks_clinic)::mldmon
      real(r4),dimension(imt,jmt,max_blocks_clinic)::ifracmon !LPF20160817
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::akmmon,aktmon,aksmon
      real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::tendmon
      real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dt_diffmon
      real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::axmon,aymon,azmon
      real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dxmon,dymon,dzmon
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::penmon
      real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dt_convmon
      !real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::hdifmon !horizontal diffusion
      real(r4),allocatable, dimension(:,:,:):: psi_euler, psi_eddy, psi
      real(r4),allocatable, dimension(:,:):: mth, mth_adv, mth_adv_iso, mth_dif
#if ( defined TIDEMIX )
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::aktidemon
      real(r4),dimension(imt,jmt,max_blocks_clinic)::wavedismon

      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::richardsonmon     !yuzp-2016/11/13
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::fztidalmon     !yuzp-2016/11/13
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp3_tidalmon     !yuzp-2016/11/13
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::ak_tide1mon     !yuzp-2016/11/19
#endif
#if ( defined CANUTOMIXOUT )
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp1_canutomon     !yuzp-2016/12/4
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp2_canutomon     !yuzp-2016/12/4
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp3_canutomon     !yuzp-2016/12/4
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp4_canutomon     !yuzp-2016/12/4
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp5_canutomon     !yuzp-2016/12/4
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp6_canutomon     !yuzp-2016/12/4
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp7_canutomon     !yuzp-2016/12/4
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp8_canutomon     !yuzp-2016/12/4

      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wk1_canutomon     !yuzp-2016/12/4
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wk2_canutomon     !yuzp-2016/12/4
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wk3_canutomon     !yuzp-2016/12/4
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wk4_canutomon     !yuzp-2016/12/4

      real(r4),dimension(imt,jmt,max_blocks_clinic)::fcor_canutomon     !yuzp-2016/12/4
      real(r4),dimension(imt,jmt,max_blocks_clinic)::fcort_canutomon     !yuzp-2016/12/4

      real(r4),dimension(imt,jmt,max_blocks_clinic)::wp10_canutomon     !yuzp-2016/12/8
      real(r4),dimension(imt,jmt,max_blocks_clinic)::wp11_canutomon     !yuzp-2016/12/8

      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp12_canutomon     !yuzp-2016/12/4
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp13_canutomon     !yuzp-2016/12/4

      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::alpha_canutomon     !yuzp-2016/12/4
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::beta_canutomon     !yuzp-2016/12/4
#endif

#if (defined ISO_TYPE_BF)
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::athkdfmon
#endif
#ifdef ISO
      real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::axmon_iso,aymon_iso,azmon_iso
      real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dxmon_iso,dymon_iso,dzmon_iso
      real(r4),dimension(imt,jmt,km,max_blocks_clinic):: wsmon_iso,vsmon_iso
#if ( defined ISOOUT )
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::vetisomon,vntisomon,vbtisomon
#endif
#endif

#if (defined SMAG_OUT)
      real(r4),dimension(imt,jmt,km)::a3mon
#endif

      real(r8) ERR_norm2mon !LPF20160823
#endif 
!for LOWREN
      real(r4), parameter :: spval =1.0e+35

end module output_mod
