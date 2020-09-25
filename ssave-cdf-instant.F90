!  CVS: $Id: ssave-cdf.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      subroutine SSAVEINS
#include <def-undef.h>
!     =================
!     output in NETcdf format
!     written by liu hai long 2001 jun
use param_mod
use pconst_mod
use dyn_mod
use tracer_mod
use forc_mod
use work_mod
use buf_mod, only:t_cpl,s_cpl,u_cpl,v_cpl,dhdx,dhdy,qheat      
use domain
use gather_scatter
use distribution
use msg_mod
use blocks
use cdf_mod
use output_netcdf !output_1dto4d
use output_mod,only: spval
!use mct_mod
!use esmf
!use seq_flds_mod
!use seq_cdata_mod
!use seq_infodata_mod
!use seq_timemgr_mod

  implicit none
#include <netcdf.inc>


      logical       :: write_restart     ! restart now
      logical:: ugrid 
      character (len=18) :: fname
      character (len=30) :: ncname
      integer :: klevel, iiday, nwmf
      integer(kind(1))   :: curr_ymd     ! Current date YYYYMMDD
      CHARACTER ( LEN =   4 ) :: ftail
      !CHARACTER ( LEN =  24 ) :: fname
      CHARACTER ( LEN =  25 ) :: fname1
      CHARACTER ( LEN =   8 ) :: dd
      CHARACTER ( LEN =   10 ) :: tt
      CHARACTER ( LEN =   5 ) :: zz
      INTEGER(r4)             :: vv(8)
      integer :: ncid1,ncid2,ncid3,ncid4,ncid5,ncid6,ncid7,ncid8,ncid9,ncid10
      real(r4),allocatable,dimension(:,:):: buffer_r4_global1 
      integer  chunks4(4),chunks1(1),chunks2(2),chunks3(3)
      integer, parameter :: ziplev=1!9
!
         number_day = iday +1 !LPF 20120816
         !number_month = month
      !if (mod(number_month,12)==0) then
         if ( number_day  > imd ) then
           if (mod(number_month,12)==0) then
            !DO LPF 20200401
            iyfm=int(number_month/12)+1
            !iyfm=iyfm+1 !20200401
            !number_month=1
            number_month=number_month+1
            !DO LPF 20200401 
            number_day=1
           else
             number_day = 1
             number_month = number_month + 1
           !if (mod(number_month,12)==0) then
           ! iyfm=iyfm+1
           end if
         end if
         nwmf= iyfm
         !DO LPF 20200401
         !mon0=number_month
         mon0=number_month-(iyfm-1)*12
         !DO LPF 20200401
         !mon0 = mod(number_month-1,12) + 1
         !mon0=number_month
         iiday = nnmonth(mon0) + iday + 1+ (nwmf-1)*365

!
     write_restart = .false.
     if ( rest_freq < 0) then
        if (number_day == 1) write_restart = .true.
     else if ( mod(iiday-1,rest_freq) == 0  ) then
        write_restart = .true.
     end if
     if(mytid==0) write(*,*)'in instant,iday=,imd=,mon0=, nwmf=, rest_freq= ',iday,imd,mon0,nwmf,rest_freq
     if(mytid==0) write(*,*)'in instant,write_restart', write_restart
!
!    if ( mod(iday,rest_freq) == 0 .or. iday == imd .or. iday == 10 .or. iday ==20) then
!    if ( mod(iday,rest_freq) == 0 .or. iday == 1 ) then
!    if ( mod(iiday,rest_freq) == 1  ) then
     if ( write_restart  ) then
!
         fname(1:8)='fort.99.'
         fname(13:13)='-'
         fname(16:16)='-'
         write(fname(14:15),'(i2.2)')mon0
         write(fname(9:12),'(i4.4)')nwmf
         write(fname(17:18),'(i2.2)') number_day
 
!        if ( number_day ==1 .and. mon0 < 12) then
!            write(fname(14:15),'(i2.2)')mon0+1
!        end if
!
!        if ( number_day ==1 .and. mon0 == 12) then
!            write(fname(14:15),'(i2.2)')mon0-11
!            write(fname(9:12),'(i4.4)')nwmf+1
!        end if

       if ( mytid == 0 ) then
!         open (17, file="rpointer.ocn", form='formatted')
!         write(17,'(a18)') fname
!         close(17)
         fname(1:8)='fort.22.' 
         write(*,*)'fname=',fname
         open(22,file=trim(out_dir)//fname,form='unformatted')
       end if

       if ( mytid == 50 ) then
         fname(1:8)='fort.23.'
         open(23,file=trim(out_dir)//fname,form='unformatted')
       end if
       if ( mytid == 20 ) then
         fname(1:8)='fort.24.'
         open(24,file=trim(out_dir)//fname,form='unformatted')
       end if
       if ( mytid == 30 ) then
         fname(1:8)='fort.25.'
         open(25,file=trim(out_dir)//fname,form='unformatted')
       end if
       if ( mytid == 40 ) then
         fname(1:8)='fort.26.'
         open(26,file=trim(out_dir)//fname,form='unformatted')
       end if
       if ( mytid == 0 ) then
         fname(1:8)='fort.27.'
         open(27,file=trim(out_dir)//fname,form='unformatted')
         fname(1:8)='fort.22.'
       end if
!================set the file name================================
      write (ftail,'(i4.4)') nwmf
      fname1(1:11)='ssh-snpshot'
      fname1(12:12)='-'
      write(fname1(13:16),'(i4.4)')nwmf
      fname1(17:17)='-'
      write(fname1(18:19),'(i2.2)')mon0
      fname1(20:20)='-'
     write(fname1(21:22),'(i2.2)') number_day
      fname1(23:25)='.nc'
!================set the file name================================

!================allocate the space================================
       allocate(buffer_r4_global1(imt_global,jmt_global), buffer_r4_local(imt,jmt,max_blocks_clinic))
       !allocate(buffer_r4_global(imt_global,jmt_global),viv_global(imt_global,jmt_global))
       allocate(buffer(imt_global,jmt_global))

!         where ( vit(:,:,1,:) > 0.5D0 )
!            buffer_r4_local = 1 
!         elsewhere
!            buffer_r4_local = spval
!         endwhere 
!         call gather_global(vit_global,buffer_r4_local, master_task,distrb_clinic) 
!mytid==0
!         where ( viv(:,:,1,:) > 0.5D0 )
!            buffer_r4_local = 1 
!         elsewhere
!            buffer_r4_local = spval
!         endwhere 
!         call gather_global(viv_global,buffer_r4_local, master_task,distrb_clinic) 
!mytid==0
         
!for ssh
       where (vit(:,:,1,:) < 0.5D0) h0 = 0.0_r8
       call gather_global(buffer, h0, master_task,distrb_clinic)
!for us
      where (viv < 0.5D0) u = 0.0_r8
      call MPI_gather(u(:,:,:,1),imt*jmt*km,MPI_PR,buffer3d,imt*jmt*km,MPI_PR,master_task,mpi_comm_ocn,ierr)
!for vs
      where (viv < 0.5D0) v = 0.0_r8
      call MPI_gather(v(:,:,:,1),imt*jmt*km,MPI_PR,buffer3d,imt*jmt*km,MPI_PR,master_task+50,mpi_comm_ocn,ierr)
!for ts
      where (vit < 0.5D0) at(:,:,:,1,:) = 0.0_r8
      call MPI_gather(at(:,:,:,1,1),imt*jmt*km,MPI_PR,buffer3d,imt*jmt*km,MPI_PR,master_task+20,mpi_comm_ocn,ierr)
!for ss
      where (vit < 0.5D0) at(:,:,:,2,:) = 0.0_r8
      call MPI_gather(at(:,:,:,2,1),imt*jmt*km,MPI_PR,buffer3d,imt*jmt*km,MPI_PR,master_task+30,mpi_comm_ocn,ierr)
!for ws
      do klevel=1,km
         where (vit(:,:,klevel,:) < 0.5D0) ws(:,:,klevel,:) = 0.0_r8
      end do
      call MPI_gather(ws(:,:,:,1),imt*jmt*km,MPI_PR,buffer3d,imt*jmt*km,MPI_PR,master_task+40,mpi_comm_ocn,ierr)

!================ssh begin================================
     if (mytid==0) then
          WRITE (22)buffer !fort22
#ifdef DEBUGXX
      iret = nf_create(fname1, NF_NETCDF4, ncid1)
      call check_err(iret)
      ugrid= .false.
      call nc_defdim(ncid1,ugrid)
       !write(*,*)'ok,define1'
      !write(*,*)'into',ncid1,iret,xadv_id,xadv_rank,xadv_dims
      call output_define(ncid1,z0_id,z0_rank,z0_dims,'ssh',time_dim,lev_dim,y_dim,x_dim)
       !write(*,*)'ok,define2'
      !write(*,*)'into3',ncid1,iret,xadv_id,xadv_rank,xadv_dims
      call output_attribute(ncid1,z0_id,z0_rank,z0_dims,'sea surface hegiht','meter','lon lat')
       !write(*,*)'ok,define3'
      iret = nf_enddef(ncid1)


         t0_cdf = 1 !need to test 
         start2(1) = 1
         start2(2) = 1
         count2(1) = lon_len
         count2(2) = lat_len
         iret = nf_put_vara_real (ncid1, lon_id,start2, count2,  lon_o)
         CALL check_err (iret)
         iret = nf_put_vara_real (ncid1, lat_id, start2, count2, lat_o)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid1, lev_id, lev)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid1, lev1_id, lev1)
         CALL check_err (iret)
         start1 (1)= 1
         count1 (1)= 1
         iret = nf_put_vara_double (ncid1, time_id,start1,count1,t0_cdf)
         CALL check_err (iret)
! store variables
         start3 (1)= 1
         start3 (2)= 1
         start3 (3)= 1
         count3 (1)= lon_len
         count3 (2)= lat_len
         count3 (3)= 1
         
         where ( abs(vit_global(:,:)) < 1.5D0 )
            buffer_r4_global1 = buffer 
         elsewhere
            buffer_r4_global1 = spval
         endwhere
!            buffer_r4_global1 = buffer 
         iret = nf_put_vara_real (ncid1,z0_id,start3, count3, buffer_r4_global1)
          CALL check_err (iret)
         iret = nf_CLOSE (ncid1)
          CALL check_err (iret)
#endif DEBUGXX

       end if !mytid==0
!================ssh and h0 Done================================
!
!================uu Done================================
      if (mytid==0) then
      
       do klevel=1,km
             call gather_global2(buffer, buffer3d(:,:,klevel,:), master_task,distrb_clinic)
             WRITE (22) buffer !writetofort22
       end do
 
       ncname(1:10)='./nc23.sh '
       ncname(11:28)=fname(1:18)
       ncname(29:30)=' &'
       call system(ncname,iret)

      end if!mytid==0
!=====================uu done===================================================


!=====================vv begin===================================================

      if (mytid==50) then

       do klevel=1,km
             call gather_global2(buffer, buffer3d(:,:,klevel,:), master_task,distrb_clinic)
             WRITE (23) buffer
       end do
       
       ncname(1:10)='./nc23.sh '
       ncname(11:28)=fname(1:18)
       ncname(29:30)=' &'
       call system(ncname,iret)

      end if
!=====================vv done===================================================

!=====================tt begin===================================================

      if (mytid==20) then

       do klevel=1,km
             call gather_global2(buffer, buffer3d(:,:,klevel,:), master_task,distrb_clinic)
             WRITE (24) buffer
       end do

       ncname(1:10)='./nc23.sh '
       ncname(11:28)=fname(1:18)
       ncname(29:30)=' &'
       call system(ncname,iret)

      end if
!=====================tt done===================================================
!=====================ss begin===================================================

      if (mytid==30) then

       do klevel=1,km
             call gather_global2(buffer, buffer3d(:,:,klevel,:), master_task,distrb_clinic)
             WRITE (25) buffer
       end do

       ncname(1:10)='./nc23.sh '
       ncname(11:28)=fname(1:18)
       ncname(29:30)=' &'
       call system(ncname,iret)

      end if
!=====================ss end===================================================
!=====================ww begin==================================================

      if (mytid==40) then

       do klevel=1,km
             call gather_global2(buffer, buffer3d(:,:,klevel,:), master_task,distrb_clinic)
             WRITE (26) buffer
       end do

       ncname(1:10)='./nc23.sh '
       ncname(11:28)=fname(1:18)
       ncname(29:30)=' &'
       call system(ncname,iret)

      end if
!=====================ww end==================================================

!=====================zonal wind stress====================================
      where (viv(:,:,1,:) < 0.5D0) su = 0.0_r8
      call gather_global(buffer, su, master_task,distrb_clinic)
!mytid.ne.0 !         where ( abs(buffer_r4_global(:,:)) < 1.5D0 )

      if (mytid==0) then
         WRITE (27)buffer

#ifdef DEBUGNC
!================taux or su begin================================
      fname1(1:11)='tx-snapshot'
      iret = nf_create(fname1, NF_NETCDF4, ncid7)
      call check_err(iret)
      ugrid= .true.
      call nc_defdim(ncid7,ugrid)
       !write(*,*)'ok,define1'
      !write(*,*)'into',ncid1,iret,xadv_id,xadv_rank,xadv_dims
      call output_define(ncid7,su_id,su_rank,su_dims,'taux',time_dim,lev_dim,y_dim,x_dim)
       !write(*,*)'ok,define2'
      !write(*,*)'into3',ncid1,iret,xadv_id,xadv_rank,xadv_dims
      call output_attribute(ncid7,su_id,su_rank,su_dims,'zonal stress','Pa','ulon ulat')
       !write(*,*)'ok,define3'
      iret = nf_enddef(ncid7)


         t0_cdf = 1 !need to test 
         start2(1) = 1
         start2(2) = 1
         count2(1) = lon_len
         count2(2) = lat_len
         iret = nf_put_vara_real (ncid7, ulon_id,start2, count2,  ulon_o)
         CALL check_err (iret)
         iret = nf_put_vara_real (ncid7, ulat_id, start2, count2, ulat_o)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid7, lev_id, lev)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid7, lev1_id, lev1)
         CALL check_err (iret)
         start1 (1)= 1
         count1 (1)= 1
         iret = nf_put_vara_double (ncid7, time_id,start1,count1,t0_cdf)
         CALL check_err (iret)
! store variables
         start3 (1)= 1
         start3 (2)= 1
         start3 (3)= 1
         count3 (1)= lon_len
         count3 (2)= lat_len
         count3 (3)= 1
         
         where ( abs(viv_global(:,:)) < 1.5D0 )
            buffer_r4_global1 = buffer 
         elsewhere
            buffer_r4_global1 = spval
         endwhere
!            buffer_r4_global1 = buffer 
         
         iret = nf_put_vara_real (ncid7,su_id,start3, count3, buffer_r4_global1)
          CALL check_err (iret)
        iret = nf_CLOSE (ncid7)
          CALL check_err (iret)
#endif DEBUGNC
      end if
!         write(*,*)'finish su'

!================taux or su begin================================
      where (viv(:,:,1,:) < 0.5D0) sv = 0.0_r8
      call gather_global(buffer, sv, master_task,distrb_clinic)
      if (mytid==0) then
         WRITE (27)buffer
       
#ifdef DEBUGNC

       fname1(1:11)='ty-snapshot'
      iret = nf_create(fname1, NF_NETCDF4, ncid8)
      call check_err(iret)
      ugrid= .true.
      call nc_defdim(ncid8,ugrid)
       !write(*,*)'ok,define1'
      !write(*,*)'into',ncid1,iret,xadv_id,xadv_rank,xadv_dims
      call output_define(ncid8,sv_id,sv_rank,sv_dims,'tauy',time_dim,lev_dim,y_dim,x_dim)
       !write(*,*)'ok,define2'
      !write(*,*)'into3',ncid1,iret,xadv_id,xadv_rank,xadv_dims
      call output_attribute(ncid8,sv_id,sv_rank,sv_dims,'meridional stress','Pa','ulon ulat')
       !write(*,*)'ok,define3'
      iret = nf_enddef(ncid8)


         t0_cdf = 1 !need to test 
         start2(1) = 1
         start2(2) = 1
         count2(1) = lon_len
         count2(2) = lat_len
         iret = nf_put_vara_real (ncid8, ulon_id,start2, count2,  ulon_o)
         CALL check_err (iret)
         iret = nf_put_vara_real (ncid8, ulat_id, start2, count2, ulat_o)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid8, lev_id, lev)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid8, lev1_id, lev1)
         CALL check_err (iret)
         start1 (1)= 1
         count1 (1)= 1
         iret = nf_put_vara_double (ncid8, time_id,start1,count1,t0_cdf)
         CALL check_err (iret)
! store variables
         start3 (1)= 1
         start3 (2)= 1
         start3 (3)= 1
         count3 (1)= lon_len
         count3 (2)= lat_len
         count3 (3)= 1
         
         where ( abs(viv_global(:,:)) < 1.5D0 )
            buffer_r4_global1 = buffer 
         elsewhere
            buffer_r4_global1 = spval
         endwhere
         
!            buffer_r4_global1 = buffer 
         iret = nf_put_vara_real (ncid8,sv_id,start3, count3, buffer_r4_global1)
          CALL check_err (iret)
        iret = nf_CLOSE (ncid8)
          CALL check_err (iret)
#endif DEBUGNC
      end if
!         write(*,*)'finish sv'

!=====================heat including sw lw sshf lthf=========
!for swv
         where (vit(:,:,1,:) < 0.5D0) swv = 0.0_r8
        call gather_global(buffer, swv, master_task,distrb_clinic)
!mytid.ne.0 !         where ( abs(buffer_r4_global(:,:)) < 1.5D0 )

     if (mytid==0) then
         WRITE (27)buffer !fort27

#ifdef DEBUGNC
       fname1(1:11)='heatsnpshot'
      iret = nf_create(fname1, NF_NETCDF4, ncid9)
      call check_err(iret)
!          write(*,*)'ok generate,ncid=',ncid, time_len
! define dimensions
         iret = nf_def_dim (ncid9, 'y', lat_len, y_dim)
         CALL check_err (iret)
!
         iret = nf_def_dim (ncid9, 'x', lon_len, x_dim)
         CALL check_err (iret)
 
         iret = nf_def_dim (ncid9, 'lev', klv, lev_dim)
         CALL check_err (iret)
 
         iret = nf_def_dim (ncid9, 'time', NF_UNLIMITED, time_dim)
         CALL check_err (iret)
       
! define variables
         lat_dims (2) = y_dim
         lat_dims (1) = x_dim
         iret = nf_def_var (ncid9, 'lat', NF_REAL, lat_rank, lat_dims, lat_id)
         CALL check_err (iret)
!
         lon_dims (2) = y_dim
         lon_dims (1) = x_dim
         iret = nf_def_var (ncid9, 'lon', NF_REAL, lon_rank, lon_dims, lon_id)
         CALL check_err (iret)
!          write(*,*)'OK define lat'
         lev_dims (1) = lev_dim
         iret = nf_def_var (ncid9, 'lev', NF_REAL, lev_rank, lev_dims, lev_id)
         CALL check_err (iret)
!
         time_dims (1) = time_dim
         iret = nf_def_var (ncid9, 'time', NF_DOUBLE, time_rank, time_dims, time_id)
         CALL check_err (iret)
!          write(*,*)'OK lev'
       
       chunks2(1) = 512!560 
       chunks2(2) = 320!356 
      iret = NF_DEF_VAR_CHUNKING(ncid9,lon_id, NF_CHUNKED, chunks2)
      CALL check_err (iret)
      iret = NF_DEF_VAR_DEFLATE (ncid9,lon_id, 1,1,ziplev)
      CALL check_err (iret)
      
      iret = NF_DEF_VAR_CHUNKING(ncid9,lat_id, NF_CHUNKED, chunks2)
      CALL check_err (iret)
      iret = NF_DEF_VAR_DEFLATE (ncid9,lat_id, 1,1,ziplev)
      CALL check_err (iret)
     !     write(*,*)'OK lat'


! define variables
         lthf_dims (3) = time_dim
         lthf_dims (2) = y_dim
         lthf_dims (1) = x_dim
         iret = nf_def_var (ncid9, 'lthf', NF_REAL, lthf_rank, lthf_dims, lthf_id)
         CALL check_err (iret)

!          write(*,*)'OK1 lthf'
         sshf_dims (3) = time_dim
         sshf_dims (2) = y_dim
         sshf_dims (1) = x_dim
         iret = nf_def_var (ncid9, 'sshf', NF_REAL, sshf_rank, sshf_dims, sshf_id)
         CALL check_err (iret)

         swv_dims (3) = time_dim
         swv_dims (2) = y_dim
         swv_dims (1) = x_dim
         iret = nf_def_var (ncid9, 'sw', NF_REAL, swv_rank, swv_dims, swv_id)
         CALL check_err (iret)

!          write(*,*)'OK1 latx'
         lwv_dims (3) = time_dim
         lwv_dims (2) = y_dim
         lwv_dims (1) = x_dim
         iret = nf_def_var (ncid9, 'lw', NF_REAL, lwv_rank, lwv_dims, lwv_id)
         CALL check_err (iret)
!          write(*,*)'OK1 latx'
         
       chunks3(1) = 1536!1786 
       chunks3(2) = 1024!1141 
       chunks3(3) = 1 
      iret = NF_DEF_VAR_CHUNKING(ncid9,lthf_id, NF_CHUNKED, chunks3)
      CALL check_err (iret)
      iret = NF_DEF_VAR_DEFLATE (ncid9,lthf_id, 1,1,ziplev)
      CALL check_err (iret)

      iret = NF_DEF_VAR_CHUNKING(ncid9,sshf_id, NF_CHUNKED, chunks3)
!      CALL check_err (iret)
      iret = NF_DEF_VAR_DEFLATE (ncid9,sshf_id, 1,1,ziplev)
      CALL check_err (iret)
      iret = NF_DEF_VAR_CHUNKING(ncid9,swv_id, NF_CHUNKED, chunks3)
      CALL check_err (iret)
      iret = NF_DEF_VAR_DEFLATE (ncid9,swv_id, 1,1,ziplev)
      CALL check_err (iret)

      iret = NF_DEF_VAR_CHUNKING(ncid9,lwv_id, NF_CHUNKED, chunks3)
      CALL check_err (iret)
      iret = NF_DEF_VAR_DEFLATE (ncid9,lwv_id, 1,1,ziplev)
      CALL check_err (iret)
! assign attributes
         iret = nf_put_att_text (ncid9, lat_id, 'long_name', 21, 'latitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid9, lat_id, 'units', 13, 'degrees_north')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid9, lon_id, 'long_name', 22, 'longitude (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid9, lon_id, 'units', 12, 'degrees_east')
         CALL check_err (iret)
         

         iret = nf_put_att_text (ncid9, lev_id, 'long_name', 18, 'depth (on T grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid9, lev_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid9, lev1_id, 'long_name', 18, 'depth (on V grids)')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid9, lev1_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid9, time_id, 'long_name', 4, 'time')
         CALL check_err (iret)
         iret = nf_put_att_text(ncid9, time_id, 'units', 21, 'days since 1001-01-01')
!         iret = nf_put_att_text (ncid, time_id, 'units', 23, 'months since 0001-01-01')
         CALL check_err (iret)
         

         iret = nf_put_att_text (ncid9, swv_id, 'long_name',9, 'shortwave')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid9, swv_id, 'units', 5, 'W/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid9, swv_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid9, swv_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)

         iret = nf_put_att_text (ncid9, lthf_id, 'long_name',6, 'latent')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid9, lthf_id, 'units', 5, 'W/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid9, lthf_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid9, lthf_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)
         
         iret = nf_put_att_text (ncid9, sshf_id, 'long_name',8, 'sensible')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid9, sshf_id, 'units', 5, 'W/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid9, sshf_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
      
         iret = nf_put_att_text (ncid9, lwv_id, 'long_name',8, 'longwave')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid9, lwv_id, 'units', 5, 'W/m^2')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid9, lwv_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid9, lwv_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)
 
!   define global attribute
         CALL date_and_time (dd,tt,zz,vv)
         iret = NF_PUT_ATT_TEXT (NCID9, NF_GLOBAL, 'title', 4, 'test')
         CALL check_err (iret)
         iret = NF_PUT_ATT_TEXT (NCID9, NF_GLOBAL, 'history', 20, tt //'  '//dd)
         CALL check_err (iret)
         iret = NF_PUT_ATT_TEXT (NCID9, NF_GLOBAL, 'source', 37, 'LASG/IAP Climate system Ocean Model 3')
         CALL check_err (iret)
! leave define mode
         iret = nf_enddef (ncid9)
         CALL check_err (iret)

!----------------------------------------------------------
!     prepare data for storing
!----------------------------------------------------------
 
         t0_cdf = 1 !need to test 
         start2(1) = 1
         start2(2) = 1
         count2(1) = lon_len
         count2(2) = lat_len
         iret = nf_put_vara_real (ncid9, lon_id,start2, count2,  lon_o)
         CALL check_err (iret)
         iret = nf_put_vara_real (ncid9, lat_id, start2, count2, lat_o)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid9, lev_id, lev)
         CALL check_err (iret)
         iret = nf_put_var_real (ncid9, lev1_id, lev1)
         CALL check_err (iret)
         start1 (1)= 1
         count1 (1)= 1
         iret = nf_put_vara_double (ncid9, time_id,start1,count1,t0_cdf)
         CALL check_err (iret)
! store variables
         start3 (1)= 1
         start3 (2)= 1
         start3 (3)= 1
         count3 (1)= lon_len
         count3 (2)= lat_len
         count3 (3)= 1
         
         where ( abs(vit_global(:,:)) < 1.5D0 )
            buffer_r4_global1 = buffer 
         elsewhere
            buffer_r4_global1 = spval
         endwhere
!            buffer_r4_global1 = buffer 
         
         iret = nf_put_vara_real (ncid9,swv_id,start3, count3, buffer_r4_global1)
          CALL check_err (iret)
#endif DEBUGNC
      end if

!for lwv
      where (vit(:,:,1,:) < 0.5D0) lwv = 0.0_r8
     call gather_global(buffer, lwv, master_task,distrb_clinic)
      
      
       if (mytid==0) then
!         WRITE (27)buffer

#ifdef DEBUGNC
         where ( abs(vit_global(:,:)) < 1.5D0 )
            buffer_r4_global1 = buffer 
         elsewhere
            buffer_r4_global1 = spval
         endwhere
!            buffer_r4_global1 = buffer 
         
         iret = nf_put_vara_real (ncid9,lwv_id,start3, count3, buffer_r4_global1)
          CALL check_err (iret)
#endif DEBUGNC
      end if
!         write(*,*)'finish lwv'
      
      where (vit(:,:,1,:) < 0.5D0) sshf = 0.0_r8
     call gather_global(buffer, sshf, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (27)buffer

#ifdef DEBUGNC
         where ( abs(vit_global(:,:)) < 1.5D0 )
            buffer_r4_global1 = buffer 
         elsewhere
            buffer_r4_global1 = spval
         endwhere
!            buffer_r4_global1 = buffer 
         
         iret = nf_put_vara_real (ncid9,sshf_id,start3, count3, buffer_r4_global1)
          CALL check_err (iret)
#endif DEBUGNC
      end if
!         write(*,*)'finish sshf'

      where (vit(:,:,1,:) < 0.5D0) lthf = 0.0_r8
     call gather_global(buffer, lthf, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (27)buffer

#ifdef DEBUGNC
         where ( abs(vit_global(:,:)) < 1.5D0 )
            buffer_r4_global1 = buffer 
         elsewhere
            buffer_r4_global1 = spval
         endwhere
!            buffer_r4_global1 = buffer 
         
         iret = nf_put_vara_real (ncid9,lthf_id,start3, count3, buffer_r4_global1)
          CALL check_err (iret)
#endif DEBUGNC
     end if
!         write(*,*)'finish lthf'
       if (mytid==0) then
#ifdef DEBUGNC
        iret = nf_CLOSE (ncid9)
          CALL check_err (iret)
#endif DEBUGNC
      end if
!============================OK heat flux==================================


      where (vit(:,:,1,:) < 0.5D0) fresh = 0.0_r8
     call gather_global(buffer, fresh, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (27)buffer
     end if
!         write(*,*)'finish fresh'
     if (mytid==0) then
         write(27) number_month, number_day
     end if
!
#ifdef COUP
      where (vit(:,:,1,:) < 0.5D0) t_cpl = 0.0_r8
     call gather_global(buffer, t_cpl, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish t_cpl'
      where (vit(:,:,1,:) < 0.5D0) s_cpl = 0.0_r8
     call gather_global(buffer, s_cpl, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish s_cpl'
      where (vit(:,:,1,:) < 0.5D0) u_cpl = 0.0_r8
     call gather_global(buffer, u_cpl, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish u_cpl'
      where (vit(:,:,1,:) < 0.5D0) v_cpl = 0.0_r8
     call gather_global(buffer, v_cpl, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish v_cpl'
      where (vit(:,:,1,:) < 0.5D0) dhdx = 0.0_r8
     call gather_global(buffer, dhdx, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish dhdx'
      where (vit(:,:,1,:) < 0.5D0) dhdy = 0.0_r8
     call gather_global(buffer, dhdy, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish dhdy'
      where (vit(:,:,1,:) < 0.5D0) qheat = 0.0_r8
     call gather_global(buffer, qheat, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish q'
#endif

     deallocate( buffer)
      if(mytid==0) then
        close(22)
        close(27)
      end if
      if(mytid==50) then
        close(23)
      end if
      if(mytid==20) then
        close(24)
      end if
      if(mytid==30) then
        close(25)
      end if
      if(mytid==40) then
        close(26)
      end if
!
  end if

      deallocate (buffer_r4_local)
      deallocate (buffer_r4_global1)
      return

      end

#if (defined NETCDF) || (defined ALL)
      SUBROUTINE check_err (iret)
#include <netcdf.inc>
      INTEGER :: iret
      IF (iret /= NF_NOERR) THEN
         PRINT *, nf_strerror (iret)
         STOP
      END IF
      END SUBROUTINE check_err
#endif

