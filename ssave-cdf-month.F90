!  CVS: $Id: ssave-cdf.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      subroutine SSAVEMON
!     =================
!     output in NETcdf format
!     written by liu hai long 2001 jun
!     =================
!     output history (netcdf) and restart (binary) files
!     remove pre-compilation NETCDF/ALL/NORMAL, remove yearly mean part 
!     keep the diagnostic part, SPMD
!     written by liu hailong & Lin Pengfei 2012 July
!     =================
#include <def-undef.h>
use precision_mod
use param_mod
use constant_mod
use pconst_mod
use output_mod
use dyn_mod
use tracer_mod
use cdf_mod
use diag_mod
use msg_mod
use domain
use distribution
use gather_scatter
      IMPLICIT none
#include <netcdf.inc>
 
#if (defined LOWRES)
      logical :: hist_output,rest_output
      CHARACTER ( LEN =   4 ) :: ftail
      CHARACTER ( LEN =  24 ) :: fname
      CHARACTER ( LEN =  15 ) :: fname1
      CHARACTER ( LEN =   8 ) :: dd
      CHARACTER ( LEN =   10 ) :: tt
      CHARACTER ( LEN =   5 ) :: zz
      INTEGER(r4)             :: vv(8)
      INTEGER :: nwmf, iblock
      real :: Num_op1 

       if(dts_accum) then
        Num_op1=float(Num_output)
        if(mytid==0) write(*,*) 'accum thermal step Num_op1',Num_op1
       else
        Num_op1=float(imd)
       endif

!
   if (diag_mth .or. diag_bsf .or. diag_msf) then 
      allocate (kmt_g(imt_global, jmt_global))
      call gather_global(kmt_g, kmt, master_task, distrb_clinic)
      if (diag_mth) then 
         allocate (mth(2,n_lat_aux_grid-1) )
         allocate (mth_adv(2,n_lat_aux_grid-1) )
         allocate (mth_adv_iso(2,n_lat_aux_grid-1) )
         allocate (mth_dif(2,n_lat_aux_grid-1) )
      end if
      if (diag_msf) then
         allocate (psi_euler(2,n_lat_aux_grid-1,km+1) )
         allocate (psi_eddy(2,n_lat_aux_grid-1,km+1) )
      end if
   end if


!---------------------------------------------------------------------
!     output monthly results
!---------------------------------------------------------------------
!    file name
 
      nwmf = iyfm
      write (ftail,'(i4.4)') nwmf
      fname1(1:5)='MMEAN'
      fname1(6:9)=ftail
      fname1(10:10)='-'
      write(fname1(11:12),'(i2.2)')mon0
      fname1(13:15)='.nc'
 
!      if (iday==imd) then
      if (mod (iyfm,io_hist)==0 ) then
         hist_output=.true.
      else
         hist_output=.false.
      endif
!
       if(mytid==0) write(*,*) 'ok inssavemon'
       
      if(dts_accum) then
        Num_op1=float(Num_output)
       else
        Num_op1=float(imd)
       endif


!--------------------------------------------------------------
!     cdf output
!--------------------------------------------------------------
!     file defination
! enter define mode
       if (mytid==0) then
         iret = nf_create (fname1, NF_CLOBBER, ncid)
         CALL check_err (iret)
          write(*,*)'ok generate,ncid=',ncid, time_len
! define dimensions
         iret = nf_def_dim (ncid, 'y', lat_len, y_dim)
         CALL check_err (iret)
!
         if (diag_msf .or. diag_mth) then
            iret = nf_def_dim (ncid, 'lat_aux', n_lat_aux_grid-1, lat_aux_dim)
            CALL check_err (iret)
            iret = nf_def_dim (ncid, 'tracer_dim', 2 , tracer_dim)
            CALL check_err (iret)
         end if
!
         iret = nf_def_dim (ncid, 'x', lon_len, x_dim)
         CALL check_err (iret)
!ZWP Define Basin
          write(*,*)'begin define basin'
         iret = nf_def_dim (ncid, 'basin', basin_len, basin_dim)
         CALL check_err (iret)
          write(*,*)'OK define basin'
!ZWP DEfine Basin 2013-10-17
 
         IF (mon0 == 12) THEN
            iret = nf_def_dim (ncid, 'lev', lev_len, lev_dim)
            CALL check_err (iret)
         ELSE
            iret = nf_def_dim (ncid, 'lev', klv, lev_dim)
            CALL check_err (iret)
         END IF
 
         iret = nf_def_dim (ncid, 'lev1', lev1_len, lev1_dim)
         CALL check_err (iret)
 
         iret = nf_def_dim (ncid, 'time', NF_UNLIMITED, time_dim)
!      iret = nf_def_dim(ncid, 'time', time_len, time_dim)
         CALL check_err (iret)
!
! define variables
         lat_dims (2) = y_dim
         lat_dims (1) = x_dim
         iret = nf_def_var (ncid, 'lat', NF_REAL, lat_rank, lat_dims, lat_id)
         CALL check_err (iret)
!
         lon_dims (2) = y_dim
         lon_dims (1) = x_dim
         iret = nf_def_var (ncid, 'lon', NF_REAL, lon_rank, lon_dims, lon_id)
         CALL check_err (iret)
         
         ulat_dims (2) = y_dim
         ulat_dims (1) = x_dim
         iret = nf_def_var (ncid, 'ulat', NF_REAL, ulat_rank, ulat_dims, ulat_id)
         CALL check_err (iret)
!
         ulon_dims (2) = y_dim
         ulon_dims (1) = x_dim
         iret = nf_def_var (ncid, 'ulon', NF_REAL, ulon_rank, ulon_dims, ulon_id)
         CALL check_err (iret)
!
          write(*,*)'OK begin lat_aux'
         if (diag_msf .or. diag_mth) then
            lat_aux_dims (1) = lat_aux_dim
            iret = nf_def_var (ncid, 'lat_aux', NF_DOUBLE, lat_aux_rank, lat_aux_dims, lat_aux_id)
         end if
          write(*,*)'OK define lat_aux'

         CALL check_err (iret)
!
!ZWP
!          write(*,*)'begin define variable basin'
!         basin_dims (1) = basin_dim
!         iret = nf_def_var (ncid, 'basin', NF_REAL, basin_rank, basin_dims, basin_id)
!         CALL check_err (iret)
!          write(*,*)'OK define variable basin'
!ZWP2013-10-17
!
         lev_dims (1) = lev_dim
         iret = nf_def_var (ncid, 'lev', NF_REAL, lev_rank, lev_dims, lev_id)
         CALL check_err (iret)
!
         lev1_dims (1) = lev1_dim
         iret = nf_def_var (ncid, 'lev1', NF_REAL, lev1_rank, lev1_dims, lev1_id)
         CALL check_err (iret)
!
         time_dims (1) = time_dim
         iret = nf_def_var (ncid, 'time', NF_DOUBLE, time_rank, time_dims, time_id)
         CALL check_err (iret)
!
         area_dims (2) = y_dim
         area_dims (1) = x_dim
         iret = nf_def_var (ncid, 'area', NF_REAL, area_rank, area_dims, area_id)
         CALL check_err (iret)

!TEST 
         z0_dims (3) = time_dim
         z0_dims (2) = y_dim
         z0_dims (1) = x_dim
         iret = nf_def_var (ncid, 'z0', NF_REAL, z0_rank, z0_dims, z0_id)
         CALL check_err (iret)
!TEST 


         !call output_define(z0_id,z0_rank,z0_dims,'ssh') !LPF


 
         ic1_dims (3) = time_dim
         ic1_dims (2) = y_dim
         ic1_dims (1) = x_dim
         iret = nf_def_var (ncid, 'ic1', NF_REAL, ic1_rank, ic1_dims, ic1_id)
         CALL check_err (iret)
 
         ic2_dims (3) = time_dim
         ic2_dims (2) = y_dim
         ic2_dims (1) = x_dim
         iret = nf_def_var (ncid, 'ic2', NF_REAL, ic2_rank, ic2_dims, ic2_id)
         CALL check_err (iret)

         net1_dims (3) = time_dim
         net1_dims (2) = y_dim
         net1_dims (1) = x_dim
         iret = nf_def_var (ncid, 'net1', NF_REAL, net1_rank, net1_dims, net1_id)
         CALL check_err (iret)
 
         net2_dims (3) = time_dim
         net2_dims (2) = y_dim
         net2_dims (1) = x_dim
         iret = nf_def_var (ncid, 'net2', NF_REAL, net2_rank, net2_dims, net2_id)
         CALL check_err (iret)

         mld_dims (3) = time_dim
         mld_dims (2) = y_dim
         mld_dims (1) = x_dim
         iret = nf_def_var (ncid, 'mld', NF_REAL, mld_rank, mld_dims, mld_id)
         CALL check_err (iret)
         
         icefrac_dims (3) = time_dim
         icefrac_dims (2) = y_dim
         icefrac_dims (1) = x_dim
         iret = nf_def_var (ncid, 'ifrac', NF_REAL, icefrac_rank, icefrac_dims, icefrac_id)
         CALL check_err (iret)

         akm_dims (4) = time_dim
         akm_dims (3) = lev_dim
         akm_dims (2) = y_dim
         akm_dims (1) = x_dim
         iret = nf_def_var (ncid, 'akm', NF_REAL, akm_rank, akm_dims, akm_id)
         CALL check_err (iret)

         akt_dims (4) = time_dim
         akt_dims (3) = lev_dim
         akt_dims (2) = y_dim
         akt_dims (1) = x_dim
         iret = nf_def_var (ncid, 'akt', NF_REAL, akt_rank, akt_dims, akt_id)
         CALL check_err (iret)

         aks_dims (4) = time_dim
         aks_dims (3) = lev_dim
         aks_dims (2) = y_dim
         aks_dims (1) = x_dim
         iret = nf_def_var (ncid, 'aks', NF_REAL, aks_rank, aks_dims, aks_id)
         CALL check_err (iret)

#if ( defined TIDEMIX )
         aktide_dims (4) = time_dim
         aktide_dims (3) = lev_dim
         aktide_dims (2) = y_dim
         aktide_dims (1) = x_dim
         iret = nf_def_var (ncid, 'aktide', NF_REAL, aktide_rank, aktide_dims, aktide_id)
         CALL check_err (iret)
#endif
#if (defined ISO_TYPE_BF)
         athkdf_dims (4) = time_dim
         athkdf_dims (3) = lev_dim
         athkdf_dims (2) = y_dim
         athkdf_dims (1) = x_dim
         iret = nf_def_var (ncid, 'athkdf', NF_REAL, athkdf_rank, athkdf_dims, athkdf_id)
         CALL check_err (iret)
#endif

#ifdef ISOOUT
          isox_dims (4) = time_dim
          isox_dims (3) = lev_dim
          isox_dims (2) = y_dim
          isox_dims (1) = x_dim
          iret = nf_def_var (ncid, 'ustar', NF_REAL, isox_rank, isox_dims,isox_id)
          CALL check_err (iret)

          isoy_dims (4) = time_dim
          isoy_dims (3) = lev_dim
          isoy_dims (2) = y_dim
          isoy_dims (1) = x_dim
          iret = nf_def_var (ncid, 'vstar', NF_REAL, isoy_rank, isoy_dims,isoy_id)
          CALL check_err (iret)

          isoz_dims (4) = time_dim
          isoz_dims (3) = lev_dim
          isoz_dims (2) = y_dim
          isoz_dims (1) = x_dim
          iret = nf_def_var (ncid, 'wstar', NF_REAL, isoz_rank, isoz_dims,isoz_id)
          CALL check_err (iret)
#endif

         if (diag_budget) then
!#ifdef BUDGETOUT
          !move another file
!#endif
        endif !diag_budget


         ts_dims (4) = time_dim
         ts_dims (3) = lev_dim
         ts_dims (2) = y_dim
         ts_dims (1) = x_dim
         iret = nf_def_var (ncid, 'ts', NF_REAL, ts_rank, ts_dims, ts_id)
         CALL check_err (iret)

         ss_dims (4) = time_dim
         ss_dims (3) = lev_dim
         ss_dims (2) = y_dim
         ss_dims (1) = x_dim
         iret = nf_def_var (ncid, 'ss', NF_REAL, ss_rank, ss_dims, ss_id)
         CALL check_err (iret)

         us_dims (4) = time_dim
         us_dims (3) = lev_dim
         us_dims (2) = y_dim
         us_dims (1) = x_dim
         iret = nf_def_var (ncid, 'us', NF_REAL, us_rank, us_dims, us_id)
         CALL check_err (iret)

         vs_dims (4) = time_dim
         vs_dims (3) = lev_dim
         vs_dims (2) = y_dim
         vs_dims (1) = x_dim
         iret = nf_def_var (ncid, 'vs', NF_REAL, vs_rank, vs_dims, vs_id)
         CALL check_err (iret)

         ws_dims (4) = time_dim
         ws_dims (3) = lev_dim
         ws_dims (2) = y_dim
         ws_dims (1) = x_dim
         iret = nf_def_var (ncid, 'ws', NF_REAL, ws_rank, ws_dims, ws_id)
         CALL check_err (iret)

         su_dims (3) = time_dim
         su_dims (2) = y_dim
         su_dims (1) = x_dim
         iret = nf_def_var (ncid, 'su', NF_REAL, su_rank, su_dims, su_id)
         CALL check_err (iret) !Uwindstress

         sv_dims (3) = time_dim
         sv_dims (2) = y_dim
         sv_dims (1) = x_dim
         iret = nf_def_var (ncid, 'sv', NF_REAL, sv_rank, sv_dims, sv_id)
         CALL check_err (iret) !Vwindstress
         lthf_dims (3) = time_dim
         lthf_dims (2) = y_dim
         lthf_dims (1) = x_dim
         iret = nf_def_var (ncid, 'lthf', NF_REAL, lthf_rank, lthf_dims, lthf_id)
         CALL check_err (iret) !latent heat flux 

         qice_dims (3) = time_dim
         qice_dims (2) = y_dim
         qice_dims (1) = x_dim
         iret = nf_def_var (ncid, 'qice', NF_REAL, qice_rank, qice_dims, qice_id)
         CALL check_err (iret) 

         fresh_dims (3) = time_dim
         fresh_dims (2) = y_dim
         fresh_dims (1) = x_dim
         iret = nf_def_var (ncid, 'fresh', NF_REAL, fresh_rank, fresh_dims, fresh_id)
         CALL check_err (iret) 

         runoff_dims (3) = time_dim
         runoff_dims (2) = y_dim
         runoff_dims (1) = x_dim
         iret = nf_def_var (ncid, 'runoff', NF_REAL, runoff_rank, runoff_dims, runoff_id)
         CALL check_err (iret)

         sshf_dims (3) = time_dim
         sshf_dims (2) = y_dim
         sshf_dims (1) = x_dim
         iret = nf_def_var (ncid, 'sshf', NF_REAL, sshf_rank, sshf_dims, sshf_id)
         CALL check_err (iret) !sensible heat flux 

         lwv_dims (3) = time_dim
         lwv_dims (2) = y_dim
         lwv_dims (1) = x_dim
         iret = nf_def_var (ncid, 'lwv', NF_REAL, lwv_rank, lwv_dims, lwv_id)
         CALL check_err (iret) !longwave 

         swv_dims (3) = time_dim
         swv_dims (2) = y_dim
         swv_dims (1) = x_dim
         iret = nf_def_var (ncid, 'swv', NF_REAL, swv_rank, swv_dims, swv_id)
         CALL check_err (iret) !shortwave 

          write(*,*)'begin psi'
         if (diag_msf) then
            psi_euler_dims (4) = time_dim
            psi_euler_dims (3) = lev1_dim
            psi_euler_dims (2) = lat_aux_dim
            psi_euler_dims (1) = basin_dim !ZWP2013-10-17
            iret = nf_def_var (ncid, 'psi_euler', NF_REAL, psi_euler_rank, psi_euler_dims, psi_euler_id)
            CALL check_err (iret)
            psi_eddy_dims (4) = time_dim
            psi_eddy_dims (3) = lev1_dim
            psi_eddy_dims (2) = lat_aux_dim
            psi_eddy_dims (1) = basin_dim !ZWP2013-10-17
            iret = nf_def_var (ncid, 'psi_eddy', NF_REAL, psi_eddy_rank, psi_eddy_dims, psi_eddy_id)
            CALL check_err (iret)
         end if
          write(*,*)'OK psi'

         if (diag_bsf) then
            bsf_dims (3) = time_dim
            bsf_dims (2) = y_dim
            bsf_dims (1) = x_dim
            iret = nf_def_var (ncid, 'bsf', NF_REAL, bsf_rank, bsf_dims, bsf_id)
            CALL check_err (iret)
         end if
 
         if (diag_mth) then
            mth_adv_dims (4) = time_dim
            mth_adv_dims (3) =  tracer_dim ! For Tracers, YYQ, Jan. 26, 2014
            mth_adv_dims (2) = lat_aux_dim
            mth_adv_dims (1) = basin_dim
            iret = nf_def_var (ncid, 'mth_adv', NF_REAL, mth_adv_rank, mth_adv_dims, mth_adv_id)
            CALL check_err (iret)
            mth_adv_iso_dims (4) = time_dim
            mth_adv_iso_dims (3) =  tracer_dim ! For Tracers, YYQ, Jan. 26, 2014
            mth_adv_iso_dims (2) = lat_aux_dim
            mth_adv_iso_dims (1) = basin_dim
            iret = nf_def_var (ncid, 'mth_adv_iso', NF_REAL, mth_adv_iso_rank, mth_adv_iso_dims, mth_adv_iso_id)
            CALL check_err (iret)
            mth_dif_dims (4) = time_dim
            mth_dif_dims (3) =  tracer_dim ! For Tracers, YYQ, Jan. 26, 2014
            mth_dif_dims (2) = lat_aux_dim
            mth_dif_dims (1) = basin_dim
            iret = nf_def_var (ncid, 'mth_dif', NF_REAL, mth_dif_rank, mth_dif_dims, mth_dif_id)
            CALL check_err (iret)
         end if

#if (defined SMAG_OUT)
         am_dims (4) = time_dim
         am_dims (3) = lev_dim
         am_dims (2) = y_dim
         am_dims (1) = x_dim
         iret = nf_def_var (ncid, 'am', NF_REAL, am_rank, am_dims, am_id)
         CALL check_err (iret)
#endif
 
! assign attributes
         iret = nf_put_att_text (ncid, lat_id, 'long_name', 21, 'latitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lat_id, 'units', 13, 'degrees_north')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'units', 12, 'degrees_east')
         CALL check_err (iret)
         
         iret = nf_put_att_text (ncid, ulat_id, 'long_name', 21, 'latitude (on U grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ulat_id, 'units', 13, 'degrees_north')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ulon_id, 'long_name', 22, 'longitude (on U grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ulon_id, 'units', 12, 'degrees_east')
         CALL check_err (iret)

!ZWP
!         iret = nf_put_att_text (ncid, basin_id, 'long_name', 11, 'Basin Index)')
!         CALL check_err (iret)
!         iret = nf_put_att_text (ncid, basin_id, 'units', 31, '1 for global and 2 for Atlantic')
!         CALL check_err (iret)
!ZWP2013-10-17

         iret = nf_put_att_text (ncid, lev_id, 'long_name', 18, 'depth (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev1_id, 'long_name', 18, 'depth (on V grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lev1_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, time_id, 'long_name', 4, 'time')
         CALL check_err (iret)
!      iret = nf_put_att_text(ncid, time_id, 'units', 21, 'days since 1001-01-01')
         iret = nf_put_att_text (ncid, time_id, 'units', 23, 'months since 0001-01-01')
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, area_id, 'long_name',16, 'T grid cell area')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, area_id, 'units', 3, 'm*m')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, area_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, area_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)

!TEST
         iret = nf_put_att_text (ncid, z0_id, 'long_name', 18, 'sea surface height')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, z0_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, z0_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, z0_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)
!TEST
!        call output_attribute(z0_id,z0_rank,z0_dims,'sea surface height ','meters','lon lat')
 
         iret = nf_put_att_text (ncid, ic1_id, 'long_name', 56, 'total of number of levels involved in convection per day')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ic1_id, 'units', 6, 'levels')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ic1_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ic1_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, ic2_id, 'long_name', 35, 'number of levels ventilated per day')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ic2_id, 'units', 6, 'levels')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ic2_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ic2_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, net1_id, 'long_name', 21, 'net surface heat flux')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, net1_id, 'units', 5, 'W/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, net1_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, net1_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, net2_id, 'long_name', 37, 'net surface salt flux include restore')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, net2_id, 'units', 8, 'kg/m^s/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, net2_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, net2_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, mld_id, 'long_name', 17, 'mixed layer depth')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, mld_id, 'units', 1, 'm')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, mld_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, mld_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)
         
         iret = nf_put_att_text (ncid, icefrac_id, 'long_name', 21, 'sea ice concentration')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, icefrac_id, 'units', 7, 'm^2/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, icefrac_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, icefrac_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, akm_id, 'long_name', 28, 'turbulent vertical viscosity')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, akm_id, 'units', 5, 'm^2/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, akm_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, akm_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, akt_id, 'long_name', 33, 'turbulent heat vertical viscosity')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, akt_id, 'units', 5, 'm^2/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, akt_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, akt_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, aks_id, 'long_name', 33, 'turbulent salt vertical viscosity')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, aks_id, 'units', 5, 'm^2/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, aks_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, aks_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)

#if ( defined TIDEMIX )
         iret = nf_put_att_text (ncid, aktide_id, 'long_name', 43, 'background vertical diffusivity due to tide')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, aktide_id, 'units', 5, 'm^2/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, aktide_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, aktide_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)
#endif
#if (defined ISO_TYPE_BF)
         iret = nf_put_att_text (ncid, athkdf_id, 'long_name', 21, 'thickness diffusivity')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, athkdf_id, 'units', 5, 'm^2/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, athkdf_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, athkdf_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)
#endif
#if (defined ISOOUT)
          iret = nf_put_att_text (ncid, isox_id, 'long_name', 26, 'eddy induced zonal current')
          CALL check_err (iret)
          iret = nf_put_att_text (ncid, isox_id, 'units', 3, 'm/s')
          CALL check_err (iret)
          iret = nf_put_att_real (ncid, isox_id, '_FillValue', NF_REAL, 1, spval)
          CALL check_err (iret)
          iret = nf_put_att_text (ncid, isox_id, 'coordinates', 9,  'ulon ulat')
          CALL check_err (iret)

          iret = nf_put_att_text (ncid, isoy_id, 'long_name', 31, 'eddy induced meridional current')
          CALL check_err (iret)
          iret = nf_put_att_text (ncid, isoy_id, 'units', 3, 'm/s')
          CALL check_err (iret)
          iret = nf_put_att_real (ncid, isoy_id, '_FillValue', NF_REAL, 1, spval)
          CALL check_err (iret)
          iret = nf_put_att_text (ncid, isoy_id, 'coordinates', 9,  'ulon ulat')
          CALL check_err (iret)
 
          iret = nf_put_att_text (ncid, isoz_id, 'long_name', 29, 'eddy induced vertical current')
          CALL check_err (iret)
          iret = nf_put_att_text (ncid, isoz_id, 'units', 3, 'm/s')
          CALL check_err (iret)
          iret = nf_put_att_real (ncid, isoz_id, '_FillValue', NF_REAL, 1, spval)
          CALL check_err (iret)
          iret = nf_put_att_text (ncid, isoz_id, 'coordinates', 9,  'ulon ulat')
          CALL check_err (iret)
#endif

         if(diag_budget) then
!move to another file
         endif !diag_budget
!
         iret = nf_put_att_text (ncid, ts_id, 'long_name', 11, 'temperature')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ts_id, 'units', 10, 'centigrade')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ts_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ts_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, ss_id, 'long_name', 8, 'salinity')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ss_id, 'units', 3, 'psu')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ss_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ss_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, us_id, 'long_name', 13, 'zonal current')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, us_id, 'units', 3, 'm/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, us_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, us_id, 'coordinates', 9,  'ulon ulat')
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, vs_id, 'long_name', 18, 'meridional current')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, vs_id, 'units', 3, 'm/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, vs_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, vs_id, 'coordinates', 9,  'ulon ulat')
         CALL check_err (iret)
 
         iret = nf_put_att_text (ncid, ws_id, 'long_name', 16, 'vertical current')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ws_id, 'units', 3, 'm/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, ws_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, ws_id, 'coordinates', 9,  'ulon ulat')
         CALL check_err (iret)
!
         iret = nf_put_att_text (ncid, su_id, 'long_name', 11, 'Uwindstress')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, su_id, 'units', 2, 'Pa')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, su_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, su_id, 'coordinates', 9,  'ulon ulat')
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, sv_id, 'long_name', 11, 'Vwindstress')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, sv_id, 'units', 2, 'Pa')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, sv_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, sv_id, 'coordinates', 9,  'ulon ulat')
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, lthf_id, 'long_name', 16, 'latent heat flux')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lthf_id, 'units', 5, 'W/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, lthf_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lthf_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)


         iret = nf_put_att_text (ncid, qice_id, 'long_name', 23, 'heat flux form freezing')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, qice_id, 'units', 5, 'W/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, qice_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, qice_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, fresh_id, 'long_name', 17, 'virtual salt flux')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, fresh_id, 'units', 8, 'kg/m^2/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, fresh_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, fresh_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, runoff_id, 'long_name', 16, 'runoff from land')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, runoff_id, 'units', 8, 'kg/m^2/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, runoff_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, runoff_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, sshf_id, 'long_name', 18, 'sensible heat flux')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, sshf_id, 'units', 5, 'W/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, sshf_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, sshf_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, lwv_id, 'long_name', 8, 'Longwave')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lwv_id, 'units', 5, 'W/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, lwv_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lwv_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)

         iret = nf_put_att_text (ncid, swv_id, 'long_name', 9, 'Shortwave')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, swv_id, 'units', 5, 'W/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, swv_id, '_FillValue', NF_REAL, 1,spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, swv_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)
!
         if (diag_msf) then
             iret = nf_put_att_text (ncid, psi_euler_id, 'long_name', 46, 'Meridioanl Stream Function for Euler Velocity')
             CALL check_err (iret)
             iret = nf_put_att_text (ncid, psi_euler_id, 'units', 8, 'Sverdrup')
             CALL check_err (iret)
             iret = nf_put_att_real (ncid, psi_euler_id, '_FillValue', NF_REAL, 1, spval)
             CALL check_err (iret)
             iret = nf_put_att_text (ncid, psi_eddy_id, 'long_name', 53, 'Meridioanl Stream Function for Eddy-Induced Velocity')
             CALL check_err (iret)
             iret = nf_put_att_text (ncid, psi_eddy_id, 'units', 8, 'Sverdrup')
             CALL check_err (iret)
             iret = nf_put_att_real (ncid, psi_eddy_id, '_FillValue', NF_REAL, 1, spval)
             CALL check_err (iret)
         end if

         if (diag_bsf) then
             iret = nf_put_att_text (ncid, bsf_id, 'long_name', 26, 'Barotropic Stream Function')
             CALL check_err (iret)
             iret = nf_put_att_text (ncid, bsf_id, 'units', 8, 'Sverdrup')
             CALL check_err (iret)
             iret = nf_put_att_real (ncid, bsf_id, '_FillValue', NF_REAL, 1, spval)
             CALL check_err (iret)
             iret = nf_put_att_text (ncid, bsf_id, 'coordinates', 9,  'ulon ulat')
             CALL check_err (iret)
         end if
!
         if (diag_mth) then
             iret = nf_put_att_text (ncid, mth_adv_id, 'long_name', 54, 'Meridional Tracer Transport (Euler Velocity Component)')
             CALL check_err (iret)
             iret = nf_put_att_text (ncid, mth_adv_id, 'units', 8, 'PW or Sv')
             CALL check_err (iret)
             iret = nf_put_att_real (ncid, mth_adv_id, '_FillValue', NF_REAL, 1, spval)
             CALL check_err (iret)
             iret = nf_put_att_text (ncid, mth_adv_iso_id, 'long_name', 61, 'Meridional Tracer Transport (Eddy-induced Velocity Component)')
             CALL check_err (iret)
             iret = nf_put_att_text (ncid, mth_adv_iso_id, 'units', 8, 'PW or Sv')
             CALL check_err (iret)
             iret = nf_put_att_real (ncid, mth_adv_iso_id, '_FillValue', NF_REAL, 1, spval)
             CALL check_err (iret)
             iret = nf_put_att_text (ncid, mth_dif_id, 'long_name', 50, 'Meridional Tracer Transport (Diffussion Component)')
             CALL check_err (iret)
             iret = nf_put_att_text (ncid, mth_dif_id, 'units', 8, 'PW or Sv')
             CALL check_err (iret)
             iret = nf_put_att_real (ncid, mth_dif_id, '_FillValue', NF_REAL, 1, spval)
             CALL check_err (iret)
         end if
!
#if (defined SMAG_OUT)
         iret = nf_put_att_text (ncid, am_id, 'long_name', 20, 'horizontal viscosity')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, am_id, 'units', 6, 'm**2/s')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, am_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, am_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)
#endif
!   define global attribute
         CALL date_and_time (dd,tt,zz,vv)
         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'title', 30, 'ocn-ice simulation by LICOM2.1')
         CALL check_err (iret)
         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'history', 20, tt //'2016Aug'//dd)
         CALL check_err (iret)
         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'source', 35, 'LASG/IAP Climate system Ocean Model')
         CALL check_err (iret)
! leave define mode
         iret = nf_enddef (ncid)
         CALL check_err (iret)

!----------------------------------------------------------
!     prepare data for storing
!----------------------------------------------------------
 
         t0_cdf = month -1 !need to test 

         start2(1) = 1
         start2(2) = 1
         count2(1) = lon_len
         count2(2) = lat_len

         if(mytid==0) then
         iret = nf_put_vara_real (ncid, lon_id,start2, count2,  lon_o)
         CALL check_err (iret)
         endif !mytid==0

         if(mytid==0) then
         iret = nf_put_vara_real (ncid, lat_id, start2, count2, lat_o)
         CALL check_err (iret)
         endif !mytid==0
         

         if(mytid==0) then
         iret = nf_put_vara_real (ncid, ulon_id,start2, count2,  ulon_o)
         CALL check_err (iret)
         endif !mytid==0

         if(mytid==0) then
         iret = nf_put_vara_real (ncid, ulat_id, start2, count2, ulat_o)
         CALL check_err (iret)
         endif !mytid==0
         
         
         if (diag_msf .or. diag_mth) then
         start1 (1)= 1
         count1 (1)= n_lat_aux_grid -1
             if(mytid==0) then
                iret = nf_put_vara_double (ncid, lat_aux_id, start1, count1, lat_aux_center)
                CALL check_err (iret)
             endif !mytid==0
         endif !mytid==0
!
         iret = nf_put_var_real (ncid, lev_id, lev)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid, lev1_id, lev1)
         CALL check_err (iret)
         start1 (1)= 1
         count1 (1)= 1
         iret = nf_put_vara_double (ncid, time_id,start1,count1,t0_cdf)
         CALL check_err (iret)
 
         endif !mytid==0

! store variables
         start3 (1)= 1
         start3 (2)= 1
         start3 (3)= 1
         count3 (1)= lon_len
         count3 (2)= lat_len
         count3 (3)= 1
 
         allocate(buffer_r4_global(imt_global,jmt_global), buffer_r4_local(imt,jmt,max_blocks_clinic))
!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = tarea
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
         if(mytid==0) then         
            iret = nf_put_vara_real (ncid,area_id,start2, count2, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0

!TEST
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = z0mon !/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
          
        if(mytid==0) then         
            iret = nf_put_vara_real (ncid,z0_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0


         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = icmon(:,:,1,:)/Num_op1 !/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
!
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,ic1_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0

!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = icmon(:,:,2,:)/Num_op1 !/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
!
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,ic2_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = netmon(:,:,1,:)/OD0CP*DZP(1) !LPF20160829 !/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
!
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,net1_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0

         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = netmon(:,:,2,:)/34.7*1.0e3*DZP(1)/OD0 !LPF20160829 !/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,net2_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = mldmon !/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,mld_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
         
          where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = ifracmon !/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,icefrac_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!3D output
         
         start4 (1)= 1
         start4 (2)= 1
         start4 (4)= 1
         count4 (1)= lon_len
         count4 (2)= lat_len
         count4 (4)= 1

!         do k = 1, klv
         do k = 1, klv-1
            start4 (3)= k
            count4 (3)= 1
!
!            where ( viv(:,:,k,:) > 0.5D0 )
            where ( viv(:,:,k+1,:) > 0.5D0 )
               buffer_r4_local = akmmon(:,:,k,:) !/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,akm_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do
!
            start4 (3)= klv
            count4 (3)= 1
               buffer_r4_local = spval
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,akm_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
!         do k = 1, klv
         do k = 1, klv-1
            start4 (3)= k
            count4 (3)= 1
!
!            where ( vit(:,:,k,:) > 0.5D0 )
            where ( vit(:,:,k+1,:) > 0.5D0 )
               buffer_r4_local = aktmon(:,:,k,:) !/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,akt_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do
!
            start4 (3)= klv
            count4 (3)= 1
               buffer_r4_local = spval
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,akt_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0

!         do k = 1, klv
         do k = 1, klv-1
            start4 (3)= k
            count4 (3)= 1
!
!            where ( vit(:,:,k,:) > 0.5D0 )
            where ( vit(:,:,k+1,:) > 0.5D0 )
               buffer_r4_local = aksmon(:,:,k,:) !/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,aks_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do

            start4 (3)= klv
            count4 (3)= 1
               buffer_r4_local = spval
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,aks_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0


#if ( defined TIDEMIX )
!         do k = 1, klv
         do k = 1, klv-1
            start4 (3)= k
            count4 (3)= 1
!
!            where ( vit(:,:,k,:) > 0.5D0 )
            where ( vit(:,:,k+1,:) > 0.5D0 )
               buffer_r4_local = aktidemon(:,:,k,:) !/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,aktide_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do

            start4 (3)= klv
            count4 (3)= 1
               buffer_r4_local = spval
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,aktide_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
#endif

#if (defined ISO_TYPE_BF)
         do k = 1, klv-1
            start4 (3)= k
            count4 (3)= 1
!
            where ( vit(:,:,k+1,:) > 0.5D0 )
               buffer_r4_local = athkdfmon(:,:,k,:) !/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,athkdf_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do

            start4 (3)= klv
            count4 (3)= 1
               buffer_r4_local = spval
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,athkdf_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
#endif


!LPF20160819
         if (diag_budget) then
!move to another file
         endif !diag_budget


!LPF20160805
#ifdef ISOOUT
!for eddy-induced v 
         do k = 1, klv
            start4 (3)= k
            count4 (3)= 1
!
            where ( vit(:,:,k,:) > 0.5D0 )
               buffer_r4_local = VNTISOMON(:,:,k,:) !/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,isoy_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do

!            start4 (3)= klv
!            count4 (3)= 1
!               buffer_r4_local = spval
!            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
!            if(mytid==0) then
!               iret = nf_put_vara_real (ncid,isoy_id,start4, count4, buffer_r4_global)
!               CALL check_err (iret)
!            endif !mytid==0

!for eddy-induced u
         do k = 1, klv
            start4 (3)= k
            count4 (3)= 1
!
            where ( vit(:,:,k,:) > 0.5D0 )
               buffer_r4_local = VETISOMON(:,:,k,:) !/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,isox_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do

!            start4 (3)= klv
!            count4 (3)= 1
!               buffer_r4_local = spval
!            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
!            if(mytid==0) then
!               iret = nf_put_vara_real (ncid,isox_id,start4, count4, buffer_r4_global)
!               CALL check_err (iret)
!            endif !mytid==0

!for eddy-induced w 
         do k = 1, klv
            start4 (3)= k
            count4 (3)= 1
!
            where ( vit(:,:,k,:) > 0.5D0 )
               buffer_r4_local = VBTISOMON(:,:,k,:) !/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,isoz_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do

!            start4 (3)= klv
!            count4 (3)= 1
!               buffer_r4_local = spval
!            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
!            if(mytid==0) then
!               iret = nf_put_vara_real (ncid,isoz_id,start4, count4, buffer_r4_global)
!               CALL check_err (iret)
!            endif !mytid==0

#endif
!LPF20160805

         do k = 1, klv
            start4 (3)= k
            count4 (3)= 1
!
            where ( vit(:,:,k,:) > 0.5D0 )
               buffer_r4_local = tsmon(:,:,k,:) !/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,ts_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do

         do k = 1, klv
            start4 (3)= k
            count4 (3)= 1
!
            where ( vit(:,:,k,:) > 0.5D0 )
               buffer_r4_local = ssmon(:,:,k,:) !/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,ss_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do

         do k = 1, klv
            start4 (3)= k
            count4 (3)= 1
!
            where ( vit(:,:,k,:) > 0.5D0 )
               buffer_r4_local = wsmon(:,:,k,:) !/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,ws_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do

         do k = 1, klv
            start4 (3)= k
            count4 (3)= 1
!
            where ( viv(:,:,k,:) > 0.5D0 )
               buffer_r4_local = usmon(:,:,k,:) !/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,us_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do

         do k = 1, klv
            start4 (3)= k
            count4 (3)= 1
!
            where ( viv(:,:,k,:) > 0.5D0 )
               buffer_r4_local = vsmon(:,:,k,:) !/float(imd)
            elsewhere
               buffer_r4_local = spval
            endwhere 
            call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
            if(mytid==0) then
               iret = nf_put_vara_real (ncid,vs_id,start4, count4, buffer_r4_global)
               CALL check_err (iret)
            endif !mytid==0
!
         end do

!Taux 

         where ( viv(:,:,1,:) > 0.5D0 )
            buffer_r4_local = sumon !/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,su_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
         where ( viv(:,:,1,:) > 0.5D0 )
            buffer_r4_local = svmon !/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,sv_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = lthfmon !/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,lthf_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = qicemon !/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,qice_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = freshmon !/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,fresh_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = runoffmon !/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,runoff_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = sshfmon !/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,sshf_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = lwvmon !/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,lwv_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = swvmon !/float(imd)
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
         if(mytid==0) then
            iret = nf_put_vara_real (ncid,swv_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
     if (diag_bsf) then
        call barosf
        if (mytid ==0 ) then
            iret = nf_put_vara_real (ncid,bsf_id,start3, count3, buffer_r4_global)
            CALL check_err (iret)
        end if
     end if
!
     if (diag_msf) then
         start4 (1)= 1
         start4 (2)= 1
         start4 (3)= 1
         start4 (4)= 1
         count4 (1)= basin_len
         count4 (2)= n_lat_aux_grid-1
         count4 (3)= lev1_len
         count4 (4)= 1 
        call msf
        if (mytid ==0 ) then
            iret = nf_put_vara_real (ncid,psi_euler_id,start4, count4, psi_euler)
            CALL check_err (iret)
            iret = nf_put_vara_real (ncid,psi_eddy_id,start4, count4, psi_eddy)
            CALL check_err (iret)
        end if
     end if
!
     if (diag_mth) then
         start4 (1)= 1
         start4 (2)= 1
         start4 (3)= 1
         start4 (4)= 1
         count4 (1)= basin_len
         count4 (2)= n_lat_aux_grid - 1
         count4 (3)= 1 
         count4 (4)= 1 
        call diag_heat_transport(1)
        if (mytid ==0 ) then
            iret = nf_put_vara_real (ncid,mth_adv_id,start4, count4, mth_adv)
            CALL check_err (iret)
            iret = nf_put_vara_real (ncid,mth_adv_iso_id,start4, count4, mth_adv_iso)
            CALL check_err (iret)
            iret = nf_put_vara_real (ncid,mth_dif_id,start4, count4, mth_dif)
            CALL check_err (iret)
        end if
         start4 (1)= 1
         start4 (2)= 1
         start4 (3)= 2
         start4 (4)= 1
         count4 (1)= basin_len
         count4 (2)= n_lat_aux_grid - 1
         count4 (3)= 1 
         count4 (4)= 1 
        call diag_heat_transport(2)
        if (mytid ==0 ) then
            iret = nf_put_vara_real (ncid,mth_adv_id,start4, count4, mth_adv)
            CALL check_err (iret)
            iret = nf_put_vara_real (ncid,mth_adv_iso_id,start4, count4, mth_adv_iso)
            CALL check_err (iret)
            iret = nf_put_vara_real (ncid,mth_dif_id,start4, count4, mth_dif)
            CALL check_err (iret)
        end if
     end if
!
     if(mytid==0) then 
         iret = nf_CLOSE (ncid)
         CALL check_err (iret)
     endif

 
         CALL mm00 
!lhl20120728      IF (mod ( (month -1),12) == 0)THEN
      IF (mod ( (month -1),io_rest) == 0)THEN
         CALL yy00
      END IF
!
      deallocate (buffer_r4_local,buffer_r4_global)
!
   if (diag_mth .or. diag_bsf .or. diag_msf) deallocate (kmt_g)
   if (diag_mth) deallocate (mth, mth_adv, mth_adv_iso, mth_dif)
   !if (diag_msf) deallocate (psi_euler,psi_eddy )
   if (diag_msf) deallocate (psi_euler )
#ifdef ISO
   if (diag_msf) deallocate (psi_eddy )
#endif

#endif
      RETURN
      END 
