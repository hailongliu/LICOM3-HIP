
module output_netcdf

!use precision_mod
!use param_mod
!use pconst_mod
!use output_mod
!use cdf_mod !, only: time_dim,y_dim,x_dim,buffer_r4_local,buffer_r4_global
!use domain
!use distribution
!use gather_scatter
!integer:: lon_dim,lat_dim,lev1_dim,lev_dim,time_dim, basin_dim,y_dim,lat_aux_dim,x_dim

      IMPLICIT NONE
!      include '/THL6/home/lhl/software/netcdf-4.1.2_ifort/include/netcdf.inc'
#include <netcdf.inc>
  public :: output1d
  public :: output_define
  public :: output_attribute
  public :: output2d
  public :: nc_defdim 
  public :: nc_defdim2d 
!  public :: output_writedata
!     error status return
      integer:: iret
      !logical:: ugrid 

!     file id
     !integer:: ncid
     character tt*8,dd*10
     integer:: lon_dim,lat_dim,lev1_dim,lev_dim,time_dim, basin_dim,y_dim,lat_aux_dim,x_dim
    real*4, parameter :: spval =1.0e+35

  !SAVE
  private

!==========================================================================
  contains
!==========================================================================
      subroutine nc_defdim(nncid,ugrid)
      !use data_licom
      use cdf_mod 
      use param_mod,only: km
      use pconst_mod,only: klv 
      integer :: nncid
      logical :: ugrid
      integer  chunks4(4),chunks1(1),chunks2(2)
      integer, parameter :: ziplev=1
!
      iret = nf_def_dim(nncid ,  'y',  lat_len,  y_dim)
      call check_err(iret)
      iret = nf_def_dim(nncid ,  'x',  lon_len,  x_dim)
      call check_err(iret)
      iret = nf_def_dim(nncid,  'lev',  km,  lev_dim)
      call check_err(iret)
      iret = nf_def_dim(nncid,  'lev1',  km+1,  lev1_dim)
      call check_err(iret)
      iret = nf_def_dim(nncid, 'time', NF_UNLIMITED, time_dim) !LPF20170727
      call check_err(iret)
       
        if(ugrid) then!
! define variables
         ulat_dims (2) = y_dim
         ulat_dims (1) = x_dim
         iret = nf_def_var (nncid, 'ulat', NF_REAL, ulat_rank, ulat_dims, ulat_id)
         CALL check_err (iret)
!
         ulon_dims (2) = y_dim
         ulon_dims (1) = x_dim
         iret = nf_def_var (nncid, 'ulon', NF_REAL, ulon_rank, ulon_dims, ulon_id)
         CALL check_err (iret)
      else!
!
       lon_dims(1) = x_dim
       lon_dims(2) = y_dim
       iret = nf_def_var(nncid,  'lon', NF_REAL,  lon_rank,  lon_dims,  lon_id)
       call check_err(iret)
       lat_dims(1) = x_dim
       lat_dims(2) = y_dim
       iret = nf_def_var(nncid,  'lat', NF_REAL,  lat_rank,  lat_dims,  lat_id)
       call check_err(iret)
      endif

      lev_dims(1) = lev_dim
      iret = nf_def_var(nncid,  'lev', NF_REAL,  lev_rank,  lev_dims,  lev_id)
      call check_err(iret)
      lev1_dims(1) = lev1_dim
      iret = nf_def_var(nncid,  'lev1', NF_REAL,  lev1_rank,  lev1_dims,  lev1_id)
      call check_err(iret)
      time_dims(1) = time_dim
      iret = nf_def_var(nncid, 'time',  NF_INT, time_rank, time_dims, time_id)
      call check_err(iret)

       chunks2(1) = 512!560
       chunks2(2) = 320!356
     if(ugrid) then!
      iret = NF_DEF_VAR_CHUNKING(nncid,ulon_id, NF_CHUNKED, chunks2)
      CALL check_err (iret)
      iret = NF_DEF_VAR_DEFLATE (nncid,ulon_id, 1,1,ziplev)
      CALL check_err (iret)
      iret = NF_DEF_VAR_CHUNKING(nncid,ulat_id, NF_CHUNKED, chunks2)
      CALL check_err (iret)
      iret = NF_DEF_VAR_DEFLATE (nncid,ulat_id, 1,1,ziplev)
      CALL check_err (iret)
     else
       iret = NF_DEF_VAR_CHUNKING(nncid,lon_id, NF_CHUNKED, chunks2)
       CALL check_err (iret)
       iret = NF_DEF_VAR_DEFLATE (nncid,lon_id, 1,1,ziplev)
       CALL check_err (iret)
      iret = NF_DEF_VAR_CHUNKING(nncid,lat_id, NF_CHUNKED, chunks2)
      CALL check_err (iret)
      iret = NF_DEF_VAR_DEFLATE (nncid,lat_id, 1,1,ziplev)
      CALL check_err (iret)
     endif

       ! chunks1(1) = 1
      !iret = NF_DEF_VAR_CHUNKING(nncid,time_id, NF_CHUNKED, chunks1)
      !CALL check_err (iret)
      !iret = NF_DEF_VAR_DEFLATE (nncid,time_id, 1,1,ziplev)
      !CALL check_err (iret)
       
      ! chunks4(1) = 1786
      ! chunks4(2) = 1141
      ! chunks4(3) = 1
      ! chunks4(4) = 1
      !iret = NF_DEF_VAR_CHUNKING(ncid,us_id, NF_CHUNKED, chunks4)
      !CALL check_err (iret)
      !iret = NF_DEF_VAR_DEFLATE (ncid,us_id, 1,1,ziplev)
      !CALL check_err (iret)


     if(ugrid) then!
! assign attributes
         iret = nf_put_att_text (nncid, ulat_id, 'long_name', 21, 'latitude (on U grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (nncid, ulat_id, 'units', 13, 'degrees_north')
         CALL check_err (iret)
         iret = nf_put_att_text (nncid, ulon_id, 'long_name', 22, 'longitude (on U grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (nncid, ulon_id, 'units', 12, 'degrees_east')
         CALL check_err (iret)
     else
       iret = nf_put_att_text(nncid, lat_id, 'long_name', 21, "latitude (on T grids)" )
       call check_err(iret)
       iret = nf_put_att_text(nncid, lat_id, 'units', 13, 'degrees_north')
       call check_err(iret)
       iret = nf_put_att_text(nncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
       call check_err(iret)
       iret = nf_put_att_text(nncid, lon_id, 'units', 12, 'degrees_east')
       call check_err(iret)
      endif
       iret = nf_put_att_text(nncid, lev_id, 'long_name', 18, 'depth (on T grids)')
       call check_err(iret)
       iret = nf_put_att_text(nncid, lev_id, 'units', 5, 'meter')
       call check_err(iret)
       iret = nf_put_att_text(nncid, lev1_id, 'long_name', 18, 'depth (on W grids)')
      call check_err(iret)
      iret = nf_put_att_text(nncid, lev1_id, 'units', 5, 'meter')
      call check_err(iret)
      iret = nf_put_att_text(nncid, time_id, 'long_name', 4, 'time')
      call check_err(iret)
      iret = nf_put_att_text(nncid, time_id, 'units', 23, 'days since 0001-01-01')
      call check_err(iret)
      iret = nf_put_att_text(nncid, time_id, 'calendar', 7, '365days')
      call check_err(iret)
!
      return
      end

      subroutine nc_defdim2d(nncid)
      !use data_licom
      use cdf_mod 
      integer :: nncid
      integer  chunks4(4),chunks1(1),chunks2(2)
      integer, parameter :: ziplev=1
!
      iret = nf_def_dim(nncid ,  'y',  lat_len,  y_dim)
      call check_err(iret)
      iret = nf_def_dim(nncid ,  'x',  lon_len,  x_dim)
      call check_err(iret)
      iret = nf_def_dim(nncid,  'lev',  lev_len,  lev_dim)
      call check_err(iret)
      iret = nf_def_dim(nncid, 'time', NF_UNLIMITED, time_dim)
      call check_err(iret)
!
      lon_dims(1) = x_dim
      lon_dims(2) = y_dim
      iret = nf_def_var(nncid,  'lon', NF_REAL,  lon_rank,  lon_dims,lon_id)
      call check_err(iret)
      lat_dims(1) = x_dim
      lat_dims(2) = y_dim
      iret = nf_def_var(nncid,  'lat', NF_REAL,  lat_rank,  lat_dims,lat_id)
      call check_err(iret)
      lev_dims(1) = lev_dim
      iret = nf_def_var(nncid,  'lev', NF_REAL,  lev_rank,  lev_dims,lev_id)
      call check_err(iret)
      time_dims(1) = time_dim
      iret = nf_def_var(nncid, 'time',  NF_INT, time_rank, time_dims,time_id)
      call check_err(iret)

       chunks2(1) = 512!560
       chunks2(2) = 320!356
      iret = NF_DEF_VAR_CHUNKING(nncid,lon_id, NF_CHUNKED, chunks2)
      CALL check_err (iret)
      iret = NF_DEF_VAR_DEFLATE (nncid,lon_id, 1,1,ziplev)
      CALL check_err (iret)
       chunks2(1) = 512!560
       chunks2(2) = 320!356
      iret = NF_DEF_VAR_CHUNKING(nncid,lat_id, NF_CHUNKED, chunks2)
      CALL check_err (iret)
      iret = NF_DEF_VAR_DEFLATE (nncid,lat_id, 1,1,ziplev)
      CALL check_err (iret)
       chunks1(1) = 1
      iret = NF_DEF_VAR_CHUNKING(nncid,time_id, NF_CHUNKED, chunks1)
      CALL check_err (iret)
      iret = NF_DEF_VAR_DEFLATE (nncid,time_id, 1,1,ziplev)
      CALL check_err (iret)

      iret = nf_put_att_text(nncid, lat_id, 'long_name', 21, "latitude (on T grids)" )
      call check_err(iret)
      iret = nf_put_att_text(nncid, lat_id, 'units', 13, 'degrees_north')
      call check_err(iret)
      iret = nf_put_att_text(nncid, lon_id, 'long_name', 22, 'longitude (on T grids)')
      call check_err(iret)
      iret = nf_put_att_text(nncid, lon_id, 'units', 12, 'degrees_east')
      call check_err(iret)
      iret = nf_put_att_text(nncid, lev_id, 'long_name', 18, 'depth (on T grids)')
      call check_err(iret)
      iret = nf_put_att_text(nncid, lev_id, 'units', 5, 'meter')
      call check_err(iret)
      iret = nf_put_att_text(nncid, time_id, 'long_name', 4, 'time')
      call check_err(iret)
      iret = nf_put_att_text(nncid, time_id, 'units', 23, 'days since 0001-01-01')
      call check_err(iret)
      iret = nf_put_att_text(nncid, time_id, 'calendar', 7, '365days')
      call check_err(iret)
!
      return
      end

!
!  CVS: $Id: output1d.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ============================================
      subroutine output1d()
!     ============================================
! output 1D data by LPF 
!      IMPLICIT NONE
!      logical :: hist_output,rest_output 
          return
          end subroutine output1d 

!     ================================================================================
     subroutine output_define(ncid,x2d_id,x2d_rank,x2d_dims,var2d_shortname,time_dim,lev1_dim,y_dim,x_dim)
!      include '/THL6/home/lhl/software/netcdf-4.1.2_ifort/include/netcdf.inc'
#include <netcdf.inc>
     !subroutine output2d(var2d,ixx,jyy,kk,lat_len,lon_len,x2d_id,ncid1,iret1,start3,count3)
!     ================================================================================
      !IMPLICIT NONE
      integer  x2d_id,ncid,iret
      integer x2d_rank 
      integer x2d_dims(x2d_rank)
      character(len=*):: var2d_shortname !,var2d_unit,var2d_coord
      integer:: lon_dim,lat_dim,lev1_dim,lev_dim,time_dim, basin_dim,y_dim,lat_aux_dim,x_dim
      integer  chunks4(4),chunks3(3)
      integer, parameter :: ziplev=1
        !write(*,*)'into',ncid,x2d_id,x2d_rank,x2d_dims,var2d_shortname
        if(x2d_rank > 3) then
         x2d_dims (4) = time_dim
         x2d_dims (3) = lev1_dim
         x2d_dims (2) = y_dim
         x2d_dims (1) = x_dim
        else
         x2d_dims (3) = time_dim
         x2d_dims (2) = y_dim
         x2d_dims (1) = x_dim
        endif
         iret =nf_def_var(ncid,var2d_shortname,NF_REAL,x2d_rank,x2d_dims,x2d_id)
          !iret =nf_def_var(ncid,trim(var2d_shortname),NF_REAL,x2d_rank,x2d_dims,x2d_id)
        !write(*,*)'into',ncid,iret,x2d_id,x2d_rank,x2d_dims,var2d_shortname
         CALL check_err (iret)
        !write(*,*)'into',ncid,iret,x2d_id,x2d_rank,x2d_dims,var2d_shortname

       chunks4(1) = 1536 !1786
       chunks4(2) = 1024 !1141
       chunks4(3) = 1
       chunks4(4) = 1
      iret = NF_DEF_VAR_CHUNKING(ncid,x2d_id, NF_CHUNKED, chunks4)
      CALL check_err (iret)
      iret = NF_DEF_VAR_DEFLATE (ncid,x2d_id, 1,1,ziplev)
      CALL check_err (iret)

             return
        end subroutine output_define 

!     ================================================================================
     subroutine output_attribute(ncid,x2d_id,x2d_rank,x2d_dims,var2d_longname,var2d_unit,var2d_coord)
!     ================================================================================
      integer :: x2d_id,ncid,iret !,ncid1,iret1
      integer:: x2d_rank,len_longname,len_unit,len_coord 
      integer:: x2d_dims(x2d_rank)
      character(len=*):: var2d_longname,var2d_unit,var2d_coord


        len_unit=len_trim(var2d_unit)
        len_longname=len_trim(var2d_longname)
        len_coord=len_trim(var2d_coord)

       !if(mytid==0) then
         iret = nf_put_att_text (ncid, x2d_id, 'long_name', len_longname, var2d_longname)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, x2d_id, 'units', len_unit, var2d_unit)
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, x2d_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, x2d_id, 'coordinates', len_coord, var2d_coord)
         CALL check_err (iret)
         
         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'title', 16, 'outputfromLICOM3')
          call check_err(iret)
         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'history', 20, tt//'  '//dd)
          call check_err(iret)
         iret = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'source', 6,&
            'LICOM3')
          call check_err(iret)

       
             return
        end subroutine output_attribute 
      

!     ================================================================================
     subroutine output2d(ncid,iret,x2d_id,x2d_rank,x2d_dims,starttmp,counttmp)
     !subroutine output2d(var2d,ixx,jyy,kk,lat_len,lon_len,x2d_id,ncid1,iret1,start3,count3)
!     ================================================================================

      integer :: x2d_id ,ncid,iret
      !integer :: kk,x2d_id,ncid1,iret1,ixx,jyy
      integer :: x2d_rank 
      !real(r4)    :: var2d(ixx,jyy,kk) !,x2d(ixx,jyy)
      integer,parameter:: ixx=3600  ! x total grid
      integer,parameter:: jyy= 2302 !y total grid
      real    :: buffer_r4_global(ixx,jyy) !,x2d(ixx,jyy)
      integer:: starttmp(x2d_rank),counttmp(x2d_rank)
      integer :: x2d_dims(x2d_rank)
      !integer:: lon_len,lat_len
          

!        if (x2d_rank >3 ) then   
!         starttmp (1)= 1
!         starttmp (2)= 1
!         starttmp (4)= 1
!         counttmp (1)= lon_len
!         counttmp (2)= lat_len
!         counttmp (4)= 1
!         starttmp (3)= k
!         counttmp (3)= 1
!        else
!         starttmp (1)= 1
!         starttmp (2)= 1
!         starttmp (3)= 1
!         counttmp (1)= lon_len
!         counttmp (2)= lat_len
!         counttmp (3)= 1
!        endif 



!         allocate(buffer_r4_global(imt_global,jmt_global), buffer_r4_local(imt,jmt,max_blocks_clinic))
!
         !where ( vit(:,:,1,:) > 0.5D0 )
         !   buffer_r4_local = var2d 
         !elsewhere
         !   buffer_r4_local = spval
         !endwhere 
         !call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
         
            iret = nf_put_vara_real (ncid,x2d_id,starttmp, counttmp, buffer_r4_global)
            CALL check_err (iret)


             return
        end subroutine output2d 



!-----------------------------------------------------------------------
!      subroutine check_err(iret)
!-----------------------------------------------------------------------
!      include '/THL6/home/lhl/software/netcdf-4.1.2_ifort/include/netcdf.inc'
!      integer iret
!      if (iret .ne. NF_NOERR) then
!      print *, nf_strerror(iret)
!      stop
!      endif
!      end subroutine check_err



end module output_netcdf
