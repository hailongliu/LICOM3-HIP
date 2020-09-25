!  CVS: $Id: inirun.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE INIRUN
!     =================
!     INITIALIZING FOR ALL PHYSICAL FIELDS
 
#include <def-undef.h>
use param_mod
use pconst_mod
use dyn_mod
use tracer_mod
use forc_mod
use work_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
use cdf_mod, only:buffer_r4_local
!use shr_msg_mod
use shr_mpi_mod
use shr_sys_mod
use buf_mod
use control_mod
use domain
use constant_mod
use gather_scatter
use POP_HaloMod
use POP_GridHorzMod
use operators
use hmix_del2
use hmix_del4
use output_mod !LPF20170831
use diag_mod

      IMPLICIT NONE

  !----- local  ------
  integer            :: fid    ! nc domain file ID
  integer            :: dimid  ! nc dimension id
  integer            :: vid    ! nc variable ID
  integer            :: rcode  ! nc return code
  integer            :: ntim   ! temporary
  integer            :: iblock   ! temporary
  integer            :: ErrorCode   ! temporary
 
#include <netcdf.inc>
!
!     Define Variables.
      integer*4,  dimension(4) :: start(4)
      integer*4,  dimension(4) :: count(4)
#if (defined HIGHRES || defined SUPHIGH)
      integer*4,  dimension(3) :: start1(3)
      integer*4,  dimension(3) :: count1(3)
#endif
      character (len=18) :: fname
      INTEGER :: NMFF
      real(r8) :: xx,yy
      integer*4   :: ncid1, iret1
#if (defined HIGHRES || defined SUPHIGH)
      integer*4 :: ncid,iret !20180928
#endif
      allocate(h0(imt,jmt,max_blocks_clinic),u(imt,jmt,km,max_blocks_clinic),v(imt,jmt,km,max_blocks_clinic), &
               at(imt,jmt,km,ntra,max_blocks_clinic))
      !allocate(atzwp(imt,jmt,km,ntra,max_blocks_clinic))


      if (mytid==0)then
          write(6,*)"BEGINNING-----INIRUN !"
      endif 
 
#ifdef COUP
  ! lihuimin, 2012.7.17, nx,ny --> imt,jmt
  allocate (t_cpl(imt,jmt,max_blocks_clinic))
  allocate (s_cpl(imt,jmt,max_blocks_clinic))
  allocate (u_cpl(imt,jmt,max_blocks_clinic))
  allocate (v_cpl(imt,jmt,max_blocks_clinic))
  allocate (dhdx(imt,jmt,max_blocks_clinic))
  allocate (dhdy(imt,jmt,max_blocks_clinic))
  allocate (Qheat   (imt,jmt,max_blocks_clinic))

  allocate (taux (imt,jmt,max_blocks_clinic))
  allocate (tauy (imt,jmt,max_blocks_clinic))
  allocate (netsw(imt,jmt,max_blocks_clinic))
  allocate (lat1 (imt,jmt,max_blocks_clinic))
  allocate (sen  (imt,jmt,max_blocks_clinic))
  allocate (lwup (imt,jmt,max_blocks_clinic))
  allocate (lwdn (imt,jmt,max_blocks_clinic))
  allocate (melth(imt,jmt,max_blocks_clinic))
  allocate (salt (imt,jmt,max_blocks_clinic))
  allocate (prec (imt,jmt,max_blocks_clinic))
  allocate (evap (imt,jmt,max_blocks_clinic))
  allocate (meltw(imt,jmt,max_blocks_clinic))
  allocate (roff (imt,jmt,max_blocks_clinic))
  allocate (iceoff (imt,jmt,max_blocks_clinic))
  allocate (snow1 (imt,jmt,max_blocks_clinic))
  allocate (ifrac(imt,jmt,max_blocks_clinic)) 
  allocate (patm (imt,jmt,max_blocks_clinic)) 
  allocate (duu10n(imt,jmt,max_blocks_clinic))

   ! if(nstart==1) then
    if(nstart>=1) then !LPF20140423 nstart==1 initilized from OBS T/S nstart=2,initialized from any ocean run
         first_step=0
    else
         first_step=1
    endif
    num_step_per_day=0
    !num_cpl=1 !6

#endif

! if (mytid==0 ) then
!    call wrap_open ('domain_licom.nc', NF_NOWRITE, fid)

!    write(6,*) 'read domain data...'
!    call shr_sys_flush(6)
!    call wrap_inq_varid(fid, 'mask', vid)
!    call wrap_get_var_int(fid,vid,mask_g)
!    do j=1,jmt_global
!    do i=1,imt_global
!       mmm(i,j) = mask_g(i,jmt_global+1-j)
!    end do
!    end do
! end if
! call scatter_global(mask_r, mmm, master_task, distrb_clinic, &
!                         field_loc_center, field_type_scalar)

     ub = 0.0D0
     vb = 0.0D0
     h0 = 0.0D0
     ubp = 0.0D0
     vbp = 0.0D0
     h0p = 0.0D0
 
     U    = 0.0D0
     V    = 0.0D0
     UP   = 0.0D0
     VP   = 0.0D0
     WS   = 0.0D0
     H0L  = 0.0D0
     H0F  = 0.0D0
     H0BL = 0.0D0
     H0BF = 0.0D0
     UTL  = 0.0D0
     UTF  = 0.0D0
     VTL  = 0.0D0
     VTF  = 0.0D0
     NET  = 0.0D0
     PXB  = 0.0D0
     PYB  = 0.0D0
     PAX  = 0.0D0
     PAY  = 0.0D0
     WHX  = 0.0D0
     WHY  = 0.0D0
     WGP  = 0.0D0
 
!     ------------------------------------------------------------------
!     Output Arrays
!     ------------------------------------------------------------------
      CALL YY00

      if(daily_accum) then
        call DAILYAVGINI
        Num_outputacc=0
      endif 
!
 
 
      IF (NSTART == 1 .or. boundary_restore > 0 ) THEN
!
      MONTH = 1
      allocate(buffer(imt_global,jmt_global))
 
     
!

#if (defined LOWRES)
!THIS IS FOR LOW
      if (mytid==0) then
       iret=nf_open('TSinitial',nf_nowrite,ncid)
       call check_err (iret)
      end if
     
!------------------------------------------------------------------
!     READ LEVITUS ANNUAL MEAN TEMPERATURE AND SALINITY
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
     
      do k=1,km !km
!
       if (mytid == 0) then
        start(1)=1 ; count(1)=imt_global
        start(2)=1 ; count(2)=jmt_global
        start(3)=k ; count(3)=1
        start(4)=1 ; count(4)=1

        iret=nf_get_vara_double(ncid,   8,start,count, buffer)
        call check_err (iret)
        do j=1 ,jmt_global/2
        do i=1, imt_global
           xx = buffer(i,j)
           yy = buffer(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx
           buffer(i,j) = yy
        end do
        end do
!
!       do j=1 ,jmt_global
!       do i=1, imt_global/2
!          xx = buffer(i,j)
!          yy = buffer(imt_global/2+i,j)
!          buffer(imt_global/2+i,j) = xx
!          buffer(i,j) = yy
!       end do
!       end do
       end if
!
      call scatter_global(at(:,:,k,1,:), buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)

!
       if (mytid == 0 ) then
        iret=nf_get_vara_double(ncid,   7,start,count, buffer)
        call check_err (iret)
        do j=1 ,jmt_global/2
        do i=1, imt_global
           xx = buffer(i,j)
           yy = buffer(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx
           buffer(i,j) = yy
        end do
        end do
!
!       do j=1 ,jmt_global
!       do i=1, imt_global/2
!          xx = buffer(i,j)
!          yy = buffer(imt_global/2+i,j)
!          buffer(imt_global/2+i,j) = xx
!          buffer(i,j) = yy
!       end do
!       end do
       end if
!
       call scatter_global(at(:,:,k,2,:), buffer,master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!
       end do !km
         iret = nf_CLOSE (ncid)
#endif


#if (defined SUPHIGH)
      ! buffer_real4
      allocate(buffer_real4(imt_global,jmt_global))
      
     if(1) then
      if (mytid==0) then
       iret=nf_open('salt_initial',nf_nowrite,ncid)
       call check_err (iret)

       iret1=nf_open('temp_initial',nf_nowrite,ncid1)
       call check_err (iret1)
      end if
      
       do k=1,km !km
!
       if (mytid == 0) then
        start1(1)=1 ; count1(1)=imt_global
        start1(2)=1 ; count1(2)=jmt_global
        start1(3)=k ; count1(3)=1

        iret=nf_get_vara_real(ncid,   4,start1,count1, buffer_real4)
        call check_err (iret)
        do j=1 ,jmt_global/2
        do i=1, imt_global
           xx = buffer_real4(i,j)
           yy = buffer_real4(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx
           buffer(i,j) = yy
        !   buffer(i,j) = yy
        end do
        end do
        !buffer(:,:)=buffer_real4(:,:)*1.0
      end if

       call scatter_global(at(:,:,k,2,:), buffer,master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
      
       if (mytid == 0) then
        !start(1)=1 ; count(1)=imt_global
        !start(2)=1 ; count(2)=jmt_global
        !start(3)=k ; count(3)=1
        !start(4)=1 ; count(4)=1
        start1(1)=1 ; count1(1)=imt_global
        start1(2)=1 ; count1(2)=jmt_global
        start1(3)=k ; count1(3)=1

        iret1=nf_get_vara_real(ncid1,   4,start1,count1, buffer_real4)
        call check_err (iret)
        do j=1 ,jmt_global/2
        do i=1, imt_global
           xx = buffer_real4(i,j)
           yy = buffer_real4(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx
           buffer(i,j) = yy
        end do
        end do
        !buffer(:,:)=buffer_real4(:,:)*1.0

      end if

       call scatter_global(at(:,:,k,1,:), buffer,master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)


     end do !KM


      if (mytid == 0) then
       iret = nf_CLOSE (ncid)
       iret1 = nf_CLOSE (ncid1)
      endif
 
     endif !1 work
!
     if(0) then
      if (mytid==0) then
       open(62,file='salt_initial',form='unformatted')
       open(63,file='temp_initial',form='unformatted')
      end if
!
!for sality 
      do k=1,km !km
       if (mytid == 0) then
         read(62) buffer_real4

       if(0) then
        do j=1 ,jmt_global/2
        do i=1, imt_global
           xx = buffer_real4(i,j)
           yy = buffer_real4(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx
           buffer(i,j) = yy
        end do
        end do
       endif !0 nowork

       !buffer= buffer_real4*VIT(:,:,K,1)1.0
       !buffer(:,:)= buffer_real4(:,:)*VIT(:,:,K,1)
       buffer= buffer_real4*1.0
      end if
       call scatter_global(at(:,:,k,2,:), buffer,master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
      at(:,:,k,2,1)=at(:,:,k,2,1)*VIT(:,:,K,1)

!for temperature
       if (mytid == 0) then
         read(63) buffer_real4
      
       if(0) then
       do j=1 ,jmt_global/2
          do i=1, imt_global
          xx = buffer_real4(i,j)
           yy = buffer_real4(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx
           buffer(i,j) = yy
        end do
        end do
       endif !0 nowork

       buffer= buffer_real4*1.0
       !buffer(:,:)= buffer_real4(:,:)*VIT(:,:,K,1)
      end if
       call scatter_global(at(:,:,k,1,:), buffer,master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)

      at(:,:,k,1,1)=at(:,:,k,1,1)*VIT(:,:,K,1)
     enddo !km

      if (mytid==0) then
       close(62)
       close(63)
      end if

     endif !0 nowork
      deallocate(buffer_real4)
      
#endif 


#if (defined HIGHRES)

      if (mytid==0) then
       iret=nf_open('salt_initial',nf_nowrite,ncid)
       call check_err (iret)

       iret1=nf_open('temp_initial',nf_nowrite,ncid1)
       call check_err (iret1)
      end if
!
      if(1)then !work-20160501
       allocate(buffer_real4(imt_global,jmt_global))
      do k=1,km !km
!
       if (mytid == 0) then
        start(1)=1 ; count(1)=imt_global
        start(2)=1 ; count(2)=jmt_global
        start(3)=k ; count(3)=1
        start(4)=1 ; count(4)=1

        iret=nf_get_vara_real(ncid,   4,start,count, buffer_real4)
        call check_err (iret)
        !do j=1 ,jmt_global/2
        !do i=1, imt_global
        !   xx = buffer(i,j)
        !   yy = buffer(i,jmt_global+1-j)
        !   buffer(i,jmt_global+1-j) = xx
        !   buffer(i,j) = yy
        !end do
        !end do
        buffer(:,:)=buffer_real4(:,:)*1000+35.0
      end if

       call scatter_global(at(:,:,k,2,:), buffer,master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
      
       if (mytid == 0) then
        start(1)=1 ; count(1)=imt_global
        start(2)=1 ; count(2)=jmt_global
        start(3)=k ; count(3)=1
        start(4)=1 ; count(4)=1

        iret1=nf_get_vara_real(ncid1,   4,start,count, buffer_real4)
        call check_err (iret)
        !do j=1 ,jmt_global/2
        !do i=1, imt_global
        !   xx = buffer(i,j)
        !   yy = buffer(i,jmt_global+1-j)
        !   buffer(i,jmt_global+1-j) = xx
        !   buffer(i,j) = yy
        !end do
        !end do

        buffer(:,:)=buffer_real4(:,:)
      end if

       call scatter_global(at(:,:,k,1,:), buffer,master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)


     end do !K
      deallocate(buffer_real4)

     endif !work-20160501

!
      if(0)then !nowork
      do k=1,km !km
!
       if (mytid == 0) then
        start(1)=1 ; count(1)=imt_global
        start(2)=1 ; count(2)=jmt_global
        start(3)=k ; count(3)=1
        start(4)=1 ; count(4)=1

        iret=nf_get_vara_double(ncid,   7,start,count, buffer)
        call check_err (iret)
        do j=1 ,jmt_global/2
        do i=1, imt_global
           xx = buffer(i,j)
           yy = buffer(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx
           buffer(i,j) = yy
        end do
        end do

      end if

       call scatter_global(at(:,:,k,2,:), buffer,master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
      
       if (mytid == 0) then
        start(1)=1 ; count(1)=imt_global
        start(2)=1 ; count(2)=jmt_global
        start(3)=k ; count(3)=1
        start(4)=1 ; count(4)=1

        iret1=nf_get_vara_double(ncid1,   7,start,count, buffer)
        call check_err (iret)
        do j=1 ,jmt_global/2
        do i=1, imt_global
           xx = buffer(i,j)
           yy = buffer(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx
           buffer(i,j) = yy
        end do
        end do

      end if

       call scatter_global(at(:,:,k,1,:), buffer,master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)


     end do

     endif !nowork
      iret = nf_CLOSE (ncid)
      iret1 = nf_CLOSE (ncid1)
#endif

 
       call POP_HaloUpdate(at(:,:,:,1,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
       call POP_HaloUpdate(at(:,:,:,2,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)





!----------------------------------------------------
!   assign 0 to land grids of TSinital
!----------------------------------------------------
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
      DO IBLOCK= 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 1,JMT
               DO I = 1,IMT
                   if(VIT(I,J,K,IBLOCK).gt.0.5) then
                   if( AT(I,J,K,1,IBLOCK).gt.50.or.AT (I,J,K,1,IBLOCK).lt.-5.0  ) &
                     write(16,*) 'AT (I,J,K,1,IBLOCK)=',i,j,k,AT (I,J,K,1,IBLOCK) 
                   if( AT(I,J,K,2,IBLOCK).gt.70.or.AT (I,J,K,2,IBLOCK).lt.0  ) &
                     write(16,*) 'AT (I,J,K,2,IBLOCK)=',i,j,k,AT (I,J,K,2,IBLOCK) 
                   !if( AT(I,J,K,1,IBLOCK).gt.50.) AT (I,J,K,1,IBLOCK) =0.0
                   !if( AT(I,J,K,1,IBLOCK).lt.-5.) AT (I,J,K,1,IBLOCK) =0.0
                   !if( AT(I,J,K,2,IBLOCK).gt.70.) AT (I,J,K,2,IBLOCK) =35.0
                   !if( AT(I,J,K,2,IBLOCK).lt.0) AT (I,J,K,2,IBLOCK) =35.0
                   endif
                  AT (I,J,K,1,IBLOCK) = AT (I,J,K,1,IBLOCK)*VIT(I,J,K,IBLOCK)
                  AT (I,J,K,2,IBLOCK) = (AT (I,J,K,2,IBLOCK)- 35.0D0)*0.001D0*VIT(I,J,K,IBLOCK)
               END DO
            END DO
         END DO
      END DO
!

!ZWP20170218
!
!      DO IBLOCK= 1, NBLOCKS_CLINIC
!         DO K = 1,KM
!            DO J = 1,JMT
!               DO I = 1,IMT
!                  ATZWP (I,J,K,1,IBLOCK) = ATZWP (I,J,K,1,IBLOCK)*VIT(I,J,K,IBLOCK)
!                  ATZWP (I,J,K,2,IBLOCK) = (ATZWP (I,J,K,2,IBLOCK)- 35.0D0)*0.001D0*VIT(I,J,K,IBLOCK)
!               END DO
!            END DO
!         END DO
!      END DO
!ZWP20170218

      RESTORE = AT
!      RESTORE = ATZWP
      deallocate(buffer)
!      deallocate(ATZWP) !ZWP20170218
!
      end if
    
      IF ( NSTART == 1) THEN 
!
         DO K = 1,KM
               ATB (:,:,K,:,:) = AT (:,:,K,:,:)
         END DO
         ATB (:,:,0,:,:) = 0.0D0
         number_day = 1
         u = 0.0d0
         v = 0.0d0
         h0= 0.0d0
!M
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)            
#ifdef COUP
            do iblock = 1, nblocks_clinic
            do j=1,jmt
            do i=1,imt
               t_cpl (i,j,iblock)  = 273.15+at(i,j,1,1,iblock)
               s_cpl (i,j,iblock)  = at(i,j,1,2,iblock)*1000.+35.
               qheat     (i,j,iblock)  = 0.0
               u_cpl (i,j,iblock)  = 0.0
               v_cpl (i,j,iblock)  = 0.0
               dhdx  (i,j,iblock)  = 0.0
               dhdy  (i,j,iblock)  = 0.0
            end do
            end do
            end do
#endif
      ELSE
 
!     ------------------------------------------------------------------
!     READ INTERMEDIATE RESULTS (fort.22/fort.21)
!     ------------------------------------------------------------------
 
       if (mytid==0) then
          open (17,file='rpointer.ocn',form='formatted')
          read(17,'(a18)') fname
          close(17)
          open(22,file=trim(out_dir)//fname,form='unformatted')
       end if
!
         allocate (buffer(imt_global,jmt_global))
!
         if (mytid==0) then
         READ (22)buffer
         end if

!        call scatter_global(h0,buffer(2:imt_global+1,:), master_task, distrb_clinic, &
         call scatter_global(h0,buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
         call POP_HaloUpdate(h0 , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
!
         do k=1,km
         if (mytid==0) then
         READ (22)buffer
         end if
!        call scatter_global(u(:,:,k,:), buffer(2:imt_global+1,:), master_task, distrb_clinic, &
         call scatter_global(u(:,:,k,:), buffer, master_task, distrb_clinic, &
                          field_loc_swcorner, field_type_vector)
         end do
         call POP_HaloUpdate(u , POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
!
         do k=1,km
         if (mytid==0) then
         READ (22)buffer
         end if
!        call scatter_global(v(:,:,k,:),buffer(2:imt_global+1,:), master_task, distrb_clinic, &
         call scatter_global(v(:,:,k,:),buffer, master_task, distrb_clinic, &
                          field_loc_swcorner, field_type_vector)
         end do
         call POP_HaloUpdate(v , POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
!
         do k=1,km
         if (mytid==0) then
         READ (22)buffer
         end if
!        call scatter_global(at(:,:,k,1,:),buffer(2:imt_global+1,:), master_task, distrb_clinic, &
         call scatter_global(at(:,:,k,1,:),buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
         end do
         call POP_HaloUpdate(at(:,:,:,1,:) , POP_haloClinic, POP_gridHorzLocCenter , &
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
!
         do k=1,km
         if (mytid==0) then
         READ (22)buffer
         end if
!        call scatter_global(at(:,:,k,2,:),buffer(2:imt_global+1,:), master_task, distrb_clinic, &
         call scatter_global(at(:,:,k,2,:),buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
         end do
         call POP_HaloUpdate(at(:,:,:,2,:) , POP_haloClinic, POP_gridHorzLocCenter , &
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
!lhl20110728 for ws
         do k=1,km
         if (mytid==0) READ (22)buffer
         end do
         if (mytid==0) READ (22)buffer !su
         if (mytid==0) READ (22)buffer !sv
         if (mytid==0) READ (22)buffer !swv
         if (mytid==0) READ (22)buffer !sshf
         if (mytid==0) READ (22)buffer !lthf
         if (mytid==0) READ (22)buffer  !fresh
        ! if (mytid==0) READ (22)buffer
        !if (mytid==0) READ (22)buffer !add lwv !noneed for new
        ! if (mytid==0) READ (22)buffer
        ! if (mytid==0) READ (22)buffer
!lhl20110728
         if (mytid == 0) then
          read(22)number_month,number_day
           if(nstart==3) then !for branch run !LPF20170621
             number_month= yearadd*12+number_month 
             number_day= dayadd+number_day
           endif
           month= number_month
         endif
      call mpi_bcast(month,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(number_month,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(number_day,1,mpi_integer,0,mpi_comm_ocn,ierr)
            write(*,*) 'number_month =',number_month,'mon0=',mon0,&
                       'number_day=',number_day,'iday=',iday
#ifdef COUP
         if (nstart==2) then
            month=(cdate/10000-1)*12+mod(cdate,10000)/100
!M
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)            
            do iblock = 1, nblocks_clinic
            do j=1,jmt
            do i=1,imt
               ! lihuimin, TODO, to be considered
               ! lihuimin, 2012.7.23, coordinate with flux_cpl, ft. yu
               t_cpl (i,j,iblock)  = 273.15+at(i,j,1,1,iblock)
               s_cpl (i,j,iblock)  = at(i,j,1,2,iblock)*1000.+35.
               ! modi end
               qheat     (i,j,iblock)  = 0.0
               u_cpl (i,j,iblock)  = 0.0
               v_cpl (i,j,iblock)  = 0.0
               dhdx  (i,j,iblock)  = 0.0
               dhdy  (i,j,iblock)  = 0.0
            end do
            end do
            end do
         else
!LPF 20120815
!for t_cpl
          if (mytid==0) then
           READ (22)buffer
          end if
!         call scatter_global(t_cpl,buffer(2:imt_global+1,:), master_task, distrb_clinic, &
          call scatter_global(t_cpl,buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!for s_cpl
          if (mytid==0) then
           READ (22)buffer
          end if
!         call scatter_global(s_cpl, buffer(2:imt_global+1,:), master_task, distrb_clinic, &
          call scatter_global(s_cpl, buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!for u_cpl
          if (mytid==0) then
           READ (22)buffer
          end if
!for v_cpl
!         call scatter_global(u_cpl, buffer(2:imt_global+1,:), master_task, distrb_clinic, &
          call scatter_global(u_cpl, buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_vector)
          if (mytid==0) then
           READ (22)buffer
          end if
!         call scatter_global(v_cpl, buffer(2:imt_global+1,:), master_task, distrb_clinic, &
          call scatter_global(v_cpl, buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_vector)

!for dhdx 
          if (mytid==0) then
           READ (22)buffer
          end if
!         call scatter_global(dhdx,buffer(2:imt_global+1,:), master_task, distrb_clinic, &
          call scatter_global(dhdx,buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!for dhdy 
          if (mytid==0) then
           READ (22)buffer
          end if
!         call scatter_global(dhdy, buffer(2:imt_global+1,:), master_task, distrb_clinic, &
          call scatter_global(dhdy, buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!for q
          if (mytid==0) then
           READ (22)buffer
          end if
!         call scatter_global(q ,buffer(2:imt_global+1,:),  master_task, distrb_clinic, &
          call scatter_global(qheat ,buffer,  master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!
         end if
!LPF 20120815
!Yu

#endif
          deallocate(buffer)
          if (mytid == 0) CLOSE(22)

       if (simple_assm.and.nstart==1) then !LPF20170821 
             number_month= yearadd*12+number_month
             number_day= dayadd+number_day
       endif
 
         NMFF = MOD (MONTH -1,12)
 
         CALL VINTEG (U,UB)
         CALL VINTEG (V,VB)
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 1,jmt
               DO I = 1,IMT
                  UP (I,J,K,IBLOCK) = U (I,J,K,IBLOCK)
                  VP (I,J,K,IBLOCK) = V (I,J,K,IBLOCK)
                  UTF (I,J,K,IBLOCK) = U (I,J,K,IBLOCK)
                  VTF (I,J,K,IBLOCK) = V (I,J,K,IBLOCK)
                  ATB (I,J,K,1,IBLOCK) = AT (I,J,K,1,IBLOCK)
                  ATB (I,J,K,2,IBLOCK) = AT (I,J,K,2,IBLOCK)
               END DO
            END DO
         END DO
       END DO
 
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = 1,jmt
            DO I = 1,IMT
               H0P (I,J,IBLOCK)= H0 (I,J,IBLOCK)
               UBP (I,J,IBLOCK)= UB (I,J,IBLOCK)
               VBP (I,J,IBLOCK)= VB (I,J,IBLOCK)
               H0F (I,J,IBLOCK)= H0 (I,J,IBLOCK)
               H0BF (I,J,IBLOCK)= H0 (I,J,IBLOCK)
            END DO
         END DO
      END DO
      END IF
!
#ifdef BIHAR
      call init_del4t
      call init_del4u
#else
      call init_del2t
      call init_del2u
#endif
!
!For get globa surface vit and viv_global
       allocate(buffer_r4_local(imt,jmt,max_blocks_clinic))
       if (mytid==0) then
        allocate(vit_global(imt_global,jmt_global))
        allocate(viv_global(imt_global,jmt_global))
       endif
         where ( vit(:,:,1,:) > 0.5D0 )
            buffer_r4_local = 1 
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(vit_global,buffer_r4_local, master_task,distrb_clinic) 
!mytid==0
         where ( viv(:,:,1,:) > 0.5D0 )
            buffer_r4_local = 1 
         elsewhere
            buffer_r4_local = spval
         endwhere 
         call gather_global(viv_global,buffer_r4_local, master_task,distrb_clinic) 
!         call gather_global(ulon_o,ulon, master_task,distrb_clinic) 
!         call gather_global(ulat_o,ulat, master_task,distrb_clinic) 
!mytid==0
          
#if (defined LOWRES)
      if (diag_mth .or. diag_msf)  then
          call init_diagnostics
      end if
#endif
!
      deallocate(buffer_r4_local)
    
      if (mytid==0)then
          write(6,*)"END-----------INIRUN !"
#ifdef COUP
          call shr_sys_flush(6)
#endif
      endif 

      RETURN

      END SUBROUTINE INIRUN
 

