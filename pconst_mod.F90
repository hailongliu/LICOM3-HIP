!  CVS: $Id: pconst_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module pconst_mod
!
#include <def-undef.h>
use precision_mod
use param_mod
!     -----------------------------------------------------------
!     Index Fields
!     -----------------------------------------------------------
!YU   real,dimension(:),allocatable:: dyr_global
      integer,dimension(jmt):: j_global
      integer,dimension(imt):: i_global
      integer :: ix,iy
      real(r8),allocatable,dimension(:,:,:,:):: vit,viv 
      real(r4),allocatable,dimension(:,:):: vit_global,viv_global 
!     real(r8),dimension(imt_global,jmt,km):: vit_1d,viv_1d
      real(r8),dimension(jmt_global):: ahv_back
!      integer,dimension(imt,jmt,max_blocks_clinic):: basin
!     integer,dimension(imt,jmt,max_blocks_clinic):: itnu
!lhl1204
      integer,dimension(imt,jmt,max_blocks_clinic):: na
      real(r8) :: dfricmx,dwndmix
      real(r8) :: Pr_number  
#ifdef BCKMEX
      real(r8) :: diff_back_eq,diff_back_coef,diff_back_coef_max  
#endif

      ! lihuimin, 2012.7.15
!     integer :: i_num   ! actual grid number in i direction in this process
!     integer :: j_num   ! actual grid number in j direction in this process
!     integer :: i_f_num ! formal grid number in i dircetion 
!     integer :: j_f_num ! formal grid number in j direction !lhl1204
!
!
!     -----------------------------------------------------------
!     Grids
!     -----------------------------------------------------------
#if (defined NETCDF) || (defined ALL)
      real(r4),dimension(imt_global):: lon
      real(r4),dimension(jmt_global):: lat
      real(r4),dimension(imt_global,jmt_global):: lon_o,lat_o
      real(r4),dimension(imt_global,jmt_global):: ulon_o,ulat_o
      real(r4),dimension(km):: lev
      real(r4),dimension(km+1):: lev1
#endif
!lhl090729
      real(r8),dimension(s_imt,s_jmt):: s_lon
      real(r8),dimension(s_imt,s_jmt):: s_lat
!lhl090729
      real(r8),dimension(km):: zkt,dzp,odzp,odzt
      real(r8),dimension(kmp1):: zkp
      real(r8),allocatable,dimension(:,:,:):: EBEA,EBEB,EBLA,EBLB,EPEA,EPEB,EPLA,EPLB
      real(r8),dimension(imt,jmt,max_blocks_clinic):: RRD1,RRD2
!lhl060506

#if ( defined SMAG)

#endif

!Yu

      real(r8),allocatable,dimension(:,:,:)::ohbt,ohbu,dzph,hbx,hby
      real(r8),allocatable,dimension(:,:,:):: SNLAT
!     real(r8),dimension(jmt):: COSU,COST
!     real(r8),dimension(:),allocatable::COSU_global,COST_global
      integer,dimension(:),allocatable:: i_start,j_start
!     real(r8),dimension(imt)::CF1,CF2,SF1,SF2
!
!
!     -----------------------------------------------------------
!     Reference T S & coefficients for calculation of d(density)
!     -----------------------------------------------------------

      REAL(r8):: TO(KM),SO(KM),C(KM,9),PO(KM)
!lhl1204

!YU
!
!
!     -----------------------------------------------------------
!     Control Parameter
!     -----------------------------------------------------------
      INTEGER:: ISOP
!
!
!     -----------------------------------------------------------
!     Phycical Parameter
!     -----------------------------------------------------------
      real(r8)::amv,ahv,ahice
!lhl1204
      real(r8),dimension(imt,jmt,km,max_blocks_clinic)::akmu,akmt
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::akt
!lhl1204
      real(r8),dimension(jmt)::AM,AH
      real(r8),dimension(imt,jmt,km,max_blocks_clinic)::am3,ah3
!lhl
      real(r8),dimension(imt,jmt,km,max_blocks_clinic)::amx,amy
!lhl
      real(r8)::gamma
#if ( defined SMAG)
      real(r8):: D0,CP,C0F,TBICE,OD0,SAG,CAG,OD0CP,ASEA, &
                    VSEA,AFB1,AFB2,AFC1,AFC2,AFT1,AFT2,KARMAN,RR
#else
      real(r8)::  D0,CP,C0F,TBICE,OD0,SAG,CAG,OD0CP,ASEA, &
                    VSEA,AFB1,AFB2,AFC1,AFC2,AFT1,AFT2
#endif
      logical :: diag_msf, diag_bsf, diag_mth, diag_budget,test_input
!
!
      CHARACTER (LEN=3):: ABMON(12)
      CHARACTER (LEN=3):: ABMON1(12)
      CHARACTER (LEN=80):: out_dir
!
      REAL(r8):: DTB,DTC,DTS,DTB2,DTC2,ONBB,ONBC,ONCC
      INTEGER:: NBB,NCC,NSS,ISB,ISC,IST,MONTH
      INTEGER:: NMONTH(12),NNMONTH(12)
      INTEGER:: number_day, number_month
!
      ! lihuimin 2012.6.18, add REFDATE
      INTEGER :: NUMBER,NSTART,IY0,IYFM,MON0,MEND,IMD,IDAY,II,JJ,IO_HIST,IO_REST,rest_freq,hist_freq,REFDATE,boundary_restore
      integer :: klv
      integer :: kvt !use the layer for simple assimilation
!
      character (len=80) :: adv_momentum, adv_tracer
      integer(kind(1))   :: curr_ymd_licom     ! Current date YYYYMMDD
      integer(kind(1))   :: curr_ymd_cpl     ! Current date YYYYMMDD
      integer(kind(1))   :: yy_licom,mm_licom,dd_licom,tod_licom     ! year, month, day
      integer(kind(1))   :: yy_cpl,mm_cpl,dd_cpl,tod_cpl     ! year, month, day
      integer(kind(1))   :: imonth,ihour,iminute,isecond,tod,ymd,ymd_sync,tod_sync     ! year, month, day
!lhl20130419
      integer:: ocn_cpl_dt,licom_cpl_dt
      integer:: yearadd,dayadd !number year and day will be added to the restart file 
      integer:: iyfmfnew !number forcing yr will be added to the restart file 
      integer:: first_step !test if first step 
      integer:: num_step_per_day !test if first step 
      integer:: Num_output,Num_outputacc !number of accumulation for output 
      integer:: Num_dailyoutput !number of accumulation for daily output 
      logical:: dts_accum !test if daily accumulation or T/S step accumula
      logical:: daily_accum !test if daily accumulation or T/S step accumula
      logical:: dailybudget_accum !test if daily accumulation or T/S step accumula
      logical:: simple_assm !using the simple assimilation(restoring) or not 

#if ( defined TIDEMIX )
       real(r8),parameter  ::local_mixing_fraction=0.33, mixing_ef=0.20, shelf_cutoff=1500., decay_scale=1000.
       real(r8),parameter  ::back_tidalmixing=1.0d-6,max_tidalmixing=1.0D-2, max_wavedis=0.1     !yuzp-2016/12/8
       real(r8),allocatable,dimension(:,:,:,:) :: fz_tide,ak_tide
       real(r8),dimension(imt,jmt,km,max_blocks_clinic) :: fztidal,richardson,wp3_tidal     !yuzp-2016/11/13
       real(r8),dimension(imt,jmt,km,max_blocks_clinic) :: ak_tide1     !yuzp-2016/11/19
#endif
#if ( defined CANUTOMIXOUT )
       real(r8),dimension(imt,jmt,km,max_blocks_clinic) ::wp1_canuto,wp2_canuto,wp3_canuto,wp4_canuto,&
                                                          wp5_canuto,wp6_canuto,wp7_canuto,wp8_canuto     !yuzp-2016/12/4

       real(r8),dimension(imt,jmt,km,max_blocks_clinic) ::wk1_canuto,wk2_canuto,wk3_canuto,wk4_canuto     !yuzp-2016/12/4

       real(r8),dimension(imt,jmt,max_blocks_clinic) ::fcor_canuto,fcort_canuto     !yuzp-2016/12/4

       real(r8),dimension(imt,jmt,km,max_blocks_clinic) ::wp12_canuto,wp13_canuto     !yuzp-2016/12/4
       real(r8),dimension(imt,jmt,max_blocks_clinic) ::wp10_canuto,wp11_canuto     !yuzp-2016/12/4
       real(r8),dimension(imt,jmt,km,max_blocks_clinic) ::alpha_canuto,beta_canuto     !yuzp-2016/12/4

#endif

end module pconst_mod
