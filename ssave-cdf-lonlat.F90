!  CVS: $Id: ssave-cdf.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      subroutine SSAVELATLON
!     =================
!     output in NETcdf format
!     =================
!     written by  Lin Pengfei 20200313
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
use grid, only: tarea
use distribution
use gather_scatter

      IMPLICIT none
#include <netcdf.inc>
 
      logical :: hist_output,rest_output
      CHARACTER ( LEN =   4 ) :: ftail
      CHARACTER ( LEN =  24 ) :: fname
      CHARACTER ( LEN =  14 ) :: fname1
      CHARACTER ( LEN =   8 ) :: dd
      CHARACTER ( LEN =   10 ) :: tt
      CHARACTER ( LEN =   5 ) :: zz
      INTEGER(r4)             :: vv(8)
      INTEGER :: nwmf, iblock
!
!---------------------------------------------------------------------
!     output monthly results
!---------------------------------------------------------------------
!    file name
 
      nwmf = iyfm
      write (ftail,'(i4.4)') nwmf
      fname1(1:11)='latlon-area'
      fname1(12:14)='.nc'
 
!      if (iday==imd) then
!
       if(mytid==0) write(*,*) 'into ssavelonlat'


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
!
          write(*,*)'OK latx'
          write(*,*)'OK define lat'

         CALL check_err (iret)
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
         
         masksurf_dims (2) = y_dim
         masksurf_dims (1) = x_dim
         iret = nf_def_var (ncid, 'masksurf', NF_REAL, masksurf_rank, masksurf_dims, masksurf_id)
         CALL check_err (iret)
 
! assign attributes
         iret = nf_put_att_text (ncid, lat_id, 'long_name', 21, 'latitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lat_id, 'units', 13, 'degrees_north')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, lon_id, 'units', 12, 'degrees_east')
         CALL check_err (iret)
         

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
!
         iret = nf_put_att_text (ncid, masksurf_id, 'long_name',16, 'T grid cell mask')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, masksurf_id, 'units', 2, 'no')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, masksurf_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, masksurf_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)
!   define global attribute
         CALL date_and_time (dd,tt,zz,vv)
         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'title', 4, 'test')
         CALL check_err (iret)
         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'history', 20, tt //'  '//dd)
         CALL check_err (iret)
         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'source', 35, 'LASG/IAP Climate system Ocean Model')
         CALL check_err (iret)
! leave define mode
         iret = nf_enddef (ncid)
         CALL check_err (iret)

!----------------------------------------------------------
!     prepare data for storing
!----------------------------------------------------------
 
         t0_cdf = 1 !need to test 

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
!
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = 1 
         elsewhere
            buffer_r4_local = 0
         endwhere 
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
         if(mytid==0) then         
            iret = nf_put_vara_real (ncid,masksurf_id,start2, count2, buffer_r4_global)
            CALL check_err (iret)
         endif !mytid==0
!
     if(mytid==0) then 
         iret = nf_CLOSE (ncid)
         CALL check_err (iret)
     endif

 
      deallocate (buffer_r4_local,buffer_r4_global)
!
      RETURN
      END 
