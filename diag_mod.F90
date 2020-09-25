module diag_mod
#include <def-undef.h>
!-------------------------------------------------------------------------------
!
! Author: Yongqiang YU  ( 1  Dec, 2003)
!
!-------------------------------------------------------------------------------
#if (defined LOWRES)
use precision_mod
use param_mod 
use constant_mod
use pconst_mod 
use output_mod 
use tracer_mod
use work_mod 
use cdf_mod 
use msg_mod
use domain
use grid
use distribution
use gather_scatter
use broadcast
!
      implicit none
!
   real (r8), dimension(:,:),public,allocatable :: ULATD_G, TLATD_G 
   integer  , dimension(:,:,:),public,allocatable :: ng_zonal
   real (r8), dimension(:),public,allocatable ::  &
      lat_aux_center,     &! cell center latitude values (degrees north)
      lat_aux_edge         ! cell edge   latitude values (degrees north)
   integer :: n_lat_aux_grid
      contains
!

!     ===================
      SUBROUTINE msf
!     ===================
!
      implicit none
!
      real(r4) :: va(imt,jmt,km,max_blocks_clinic), vvv (imt,jmt,km, max_blocks_clinic)
      integer :: nin,nta(jmt_global), iblock

!      integer :: test_flag
!
      allocate (work_1(imt,jmt,km,max_blocks_clinic),psi(2,jmt_global,km+1))!ZWP2013-10-17,work_2(imt,jmt,km,max_blocks_clinic))
!
!
      psi = 0.0_r4
      psi_euler = 0.0_r4 !zwp add 2013-10-15
      psi_eddy = 0.0_r4 !zwp add 2013-10-15
      buffer_r4_local = 0.0_r4 !zwp add 2013-10-15
      buffer_r4_global = 0.0_r4 !zwp add 2013-10-15
!
      DO  NIN=1,2
!
      IF( NIN == 1 )then
!
!$OMP PARALLEL DO PRIVATE (iblock,K,J,I)
      do iblock = 1, nblocks_clinic
      do k=1,km
      do j=1,jmt
      do i=1,imt
         work_1(i,j,k,iblock)=vit(i,j,k,iblock)
      enddo
      enddo
      enddo
      enddo
!
      ELSE
!$OMP PARALLEL DO PRIVATE (iblock,K,J,I)
      do iblock = 1, nblocks_clinic
      do k=1,km
      do j=1,jmt
      do i=1,imt
         if (basin(i,j,iblock)==1.or.basin(i,j,iblock)==2) then
            work_1(i,j,k,iblock)=vit(i,j,k,iblock)
         else
            work_1(i,j,k,iblock)=0.0
         end if
      enddo
      enddo
      enddo
      enddo
!
!      write (*, *) "test_flag = ", test_flag
      ENDIF
!
!zwp2013-10-15      enddo
!
!
      if(NIN==1)then 
      DO K=2,km
         buffer_r4_local = work_1(:,:,k,:)*wsmon(:,:,k,:)*tarea(:,:,:) !/float(imd)
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
         if (mytid == 0 ) then
            DO jj=n_lat_aux_grid-2, 1, -1
            DO j=1,jmt_global
            DO I=1,imt_global
               if (tlatd_g(i,j) < lat_aux_edge(jj) .and. tlatd_g(i,j) >= lat_aux_edge(jj+1) ) then
                  psi_euler(NIN,JJ,K)=psi_euler(NIN,JJ,k)+buffer_r4_global(I,J)*1.0E-6
               end if
            ENDDO
            ENDDO
               psi_euler(NIN,JJ,K)=psi_euler(NIN,JJ+1,k)+ psi_euler(NIN,JJ,k)
            ENDDO
         end if
      ENDDO
#ifdef ISO
      DO K=2,km
         buffer_r4_local = work_1(:,:,k,:)*wsmon_iso(:,:,k,:)*tarea(:,:,:) !/float(imd)
         call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic)
         if (mytid == 0 ) then
            DO jj=n_lat_aux_grid-2, 1, -1
            DO j=1,jmt_global
            DO I=1,imt_global
               if (tlatd_g(i,j) < lat_aux_edge(jj) .and. tlatd_g(i,j) >= lat_aux_edge(jj+1) ) then
                  psi_eddy(NIN,JJ,K)=psi_eddy(NIN,JJ,k)+buffer_r4_global(I,J)*1.0E-6
               end if
            ENDDO
            ENDDO
               psi_eddy(NIN,JJ,K)=psi_eddy(NIN,JJ+1,k)+ psi_eddy(NIN,JJ,k)
            ENDDO
         end if
      END DO
#endif

      else

      do k=2,km+1
         buffer_r4_local    = work_1(:,:,k-1,:)*vsmon(:,:,k-1,:)*dzp(k-1)*dxu !/float(imd)
!         if ( mytid > 15 .and. mytid < 23) then
!         write(120+mytid,*) "K=",K
!         write(120+mytid,*) (buffer_r4_local(i,22,1),i=3,imt-2)
!         write(160+mytid,*) "K=",K
!         write(160+mytid,*) (work_1(i,22,k-1,1),i=3,imt-2)
!         write(200+mytid,*) "K=",K
!         write(200+mytid,*) (vsmon(i,22,k-1,1),i=3,imt-2)
!         end if
         call gather_global(buffer_r4_global,buffer_r4_local,master_task,distrb_clinic)
!         if (mytid ==0 ) then
!         write(240+mytid,*) "K=",K
!         write(240+mytid,*) (buffer_r4_global(i,130),i=1,imt_global)
!         end if

         if(mytid == 0)then
            do i=1,imt_global
               psi(NIN,130,k)=psi(NIN,130,k)+buffer_r4_global(i,130)*1.0E-6
            enddo   
               psi(NIN,130,k) = psi(NIN,130,k-1)+psi(NIN,130,k)
               psi_euler(NIN,n_lat_aux_grid-88,k) = psi(NIN,130,k)
               write(*,*) "psi =",psi(NIN,130,k)
         end if 
      enddo

!      close(120+mytid)
!      close(160+mytid)
!      close(200+mytid)
       

      do k=2,km
         buffer_r4_local = work_1(:,:,k,:)*wsmon(:,:,k,:)*tarea(:,:,:) !/float(imd)
         call gather_global(buffer_r4_global,buffer_r4_local,master_task,distrb_clinic)

      if(mytid == 0)then
      do jj = n_lat_aux_grid-89,1,-1
      do j=1,jmt_global
      do i=1,imt_global
         if(tlatd_g(i,j) < lat_aux_edge(jj) .and. tlatd_g(i,j) >= lat_aux_edge(jj+1) ) then
            psi_euler(NIN,jj,k)=psi_euler(NIN,jj,k)+buffer_r4_global(i,j)*1.0E-6
         end if
      enddo
      enddo
            psi_euler(NIN,jj,k)=psi_euler(NIN,jj+1,k)+ psi_euler(NIN,jj,k)
      enddo
      end if
      enddo

#ifdef ISO
      psi = 0.0_r4
 
      do k=2,km+1
         buffer_r4_local    = work_1(:,:,k-1,:)*vsmon_iso(:,:,k-1,:)*dzp(k-1)*dxu !/float(imd)
         call gather_global(buffer_r4_global,buffer_r4_local,master_task,distrb_clinic)
         
         if(mytid == 0)then
            do i=1,imt_global
               psi(NIN,130,k)=psi(NIN,130,k)+buffer_r4_global(i,130)*1.0E-6
            enddo
               psi(NIN,130,k) = psi(NIN,130,k-1)+psi(NIN,130,k)
               psi_eddy(NIN,n_lat_aux_grid-88,k) = psi(NIN,130,k)
               write(*,*)"psi_eddy=",psi_eddy(NIN,n_lat_aux_grid-88,k)
         end if 
       enddo

       do k=2,km
          buffer_r4_local = work_1(:,:,k,:)*wsmon_iso(:,:,k,:)*tarea(:,:,:) !/float(imd)
          call gather_global(buffer_r4_global,buffer_r4_local,master_task,distrb_clinic)

          if(mytid == 0)then
             do jj = n_lat_aux_grid-89,1,-1
             do j=1,jmt_global
             do i=1,imt_global
             if(tlatd_g(i,j) < lat_aux_edge(jj) .and. tlatd_g(i,j) >= lat_aux_edge(jj+1) ) then
                psi_eddy(NIN,jj,k)=psi_eddy(NIN,jj,k)+buffer_r4_global(i,j)*1.0E-6
             end if
             enddo
             enddo
                psi_eddy(NIN,jj,k)=psi_eddy(NIN,jj+1,k)+ psi_eddy(NIN,jj,k)
             enddo
         end if
      enddo 
#endif
      END IF
 
!
!     if (mytid == 0)  then
!       where ( abs(psi(NIN,:,:)) < 1.0e-25) psi(NIN,:,:) = spval
!     end if
!
     if (mytid == 0.and.NIN==2)  then
      do jj =1, n_lat_aux_grid-1
      !do jj = n_lat_aux_grid-89,1,-1
        if(lat_aux_edge(jj+1)<=-31.5 ) then
         do k=1,km
          psi_eddy(NIN,jj,k)=spval
          psi_euler(NIN,jj,k)=spval 
         enddo  
        endif
      enddo
     endif
 
      enddo !end NIN zwp add 2013-10-15

      where ( ng_zonal == 0 )  
        psi_eddy= spval
        psi_euler= spval
      end where
!
      deallocate (work_1,psi)!ZWP 2013-10-172,work_2)
!
      RETURN
      END SUBROUTINE MSF

!
!====================================
      SUBROUTINE BAROSF
!====================================
      implicit none
!
      integer :: iblock
!
      buffer_r4_local= 0.0_r4
      buffer_r4_global= 0.0_r4
!
      do iblock = 1, nblocks_clinic
      do k=1,km
      do j=2,jmt-1
      do i=1,imt
        buffer_r4_local(i,j,iblock)=buffer_r4_local(i,j,iblock)+viv(i,j,k,iblock)*usmon(i,j,k,iblock)* &
                                    1.0e-6*dzp(k)*dyu(i,j,iblock) !/float(imd)
      enddo
      enddo
      enddo
      enddo
!
      call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
!
      if ( mytid == 0 ) then
         do j=2,jmt_global
         do i=1,imt_global
            buffer_r4_global(i,j)= buffer_r4_global(i,j-1) + buffer_r4_global(i,j)
         end do
         end do
      end if
!
      where ( kmt_g < 1 )  buffer_r4_global = spval
!
      return
      END SUBROUTINE BAROSF
!

!========================================
      SUBROUTINE diag_heat_transport(NNN)
!========================================
!
      implicit none
      integer :: NNN,NIN,iblock
      real(r4) ,dimension(imt,jmt,km,max_blocks_clinic) :: ttpp
      real(r4) ,dimension(2,jmt_global)::tmp 
      allocate (work1_g(imt_global,jmt_global), work_1(imt,jmt,km,max_blocks_clinic))
! 
      buffer_r4_global = 0.0_r4
 
      mth_adv = 0.0_r4
      mth_dif = 0.0_r4
#ifdef ISO
      mth_adv_iso = 0.0_r4
#endif
!

! vertical averaged

     do nin = 1, 2
!
      IF( NIN == 1 )then
!
!$OMP PARALLEL DO PRIVATE (iblock,K,J,I)
      do iblock = 1, nblocks_clinic
      do k=1,km
      do j=1,jmt
      do i=1,imt
         work_1(i,j,k,iblock)=vit(i,j,k,iblock)
      enddo
      enddo
      enddo
      enddo
!
      ELSE
!$OMP PARALLEL DO PRIVATE (iblock,K,J,I)
      do iblock = 1, nblocks_clinic
      do k=1,km
      do j=1,jmt
      do i=1,imt
         if (basin(i,j,iblock)==1.or.basin(i,j,iblock)==2) then
            work_1(i,j,k,iblock)=vit(i,j,k,iblock)
         else
            work_1(i,j,k,iblock)=0.0
         end if
      enddo
      enddo
      enddo
      enddo
!
      ENDIF
!
     IF(NIN==1)then

     ttpp=0.0_r4
     work1_g = 0.0_r4
!$OMP PARALLEL DO PRIVATE (iblock,K,J,I)
      do iblock= 1, nblocks_clinic
      do k=1,km
      do j=3,jmt-2
      do i=3,imt-2
          ttpp(i,j,k,iblock)=(axmon(i,j,k,nnn,iblock)+aymon(i,j,k,nnn,iblock)+azmon(i,j,k,nnn,iblock)) & 
                             *work_1(i,j,k,iblock)*dzp(k)*tarea(i,j,iblock)
      enddo
      enddo
      enddo
      enddo
!
      do k=1,km
         call gather_global(buffer_r4_global,ttpp(:,:,k,:), master_task,distrb_clinic) 
         if (mytid == 0 ) then
         do j=1,jmt_global
         do i=1,imt_global
            work1_g(i,j)=work1_g(i,j)+ buffer_r4_global(i,j)
         enddo
         enddo
         end if
      enddo
   
!
    if (mytid == 0 ) then
            DO jj=n_lat_aux_grid-2, 1, -1
            DO j=1,jmt_global
            DO I=1,imt_global
               if (tlatd_g(i,j) < lat_aux_edge(jj) .and. tlatd_g(i,j) >= lat_aux_edge(jj+1) ) then
                  mth_adv(NIN,JJ)=mth_adv(NIN,JJ)-work1_g(I,J)
               end if
            ENDDO
            ENDDO
               mth_adv(NIN,JJ)=mth_adv(NIN,JJ+1)+ mth_adv(NIN,JJ)
            ENDDO
    end if
!
    ELSE
    tmp = 0.0_r4
    ttpp = 0.0_r4
    work1_g = 0.0_r4

!$OMP PARALLEL DO PRIVATE (iblock,K,J,I)
      do iblock = 1, nblocks_clinic
      do k =1,km
      do j= 2, jmt-1
      do i= 2, imt-1
         if(basin(i,j,iblock)==1.or.basin(i,j,iblock)==2)then
         if(NNN==1)then
            ttpp(i,j,k,iblock)=vsmon(i,j,k,iblock)*(-1)*viv(i,j,k,iblock)*0.25*&
                      (tsmon(i  ,j+1,k,iblock)+tsmon(i  ,j,k,iblock)+&
                       tsmon(i-1,j+1,k,iblock)+tsmon(i-1,j,k,iblock))*dzp(k)*dxu(i,j,iblock)
         else
            ttpp(i,j,k,iblock)=vsmon(i,j,k,iblock)*(-1)*viv(i,j,k,iblock)*0.25*0.001D0*   &
                      ((ssmon(i  ,j+1,k,iblock)+ssmon(i  ,j,k,iblock)+&
                       ssmon(i-1,j+1,k,iblock)+ssmon(i-1,j,k,iblock))-140.0D0)*dzp(k)*dxu(i,j,iblock)
         end if
         end if
      enddo
      enddo
      enddo
      enddo

      do k=1,km
         call gather_global(buffer_r4_global,ttpp(:,:,k,:), master_task,distrb_clinic)
         if (mytid == 0 ) then
             do i=1,imt_global
                work1_g(i,130)=work1_g(i,130)+ buffer_r4_global(i,130)
             enddo
         end if
      enddo

      if (mytid == 0 ) then
          do i=1,imt_global
             tmp(NIN,130)= tmp(NIN,130)-work1_g(i,130) 
          enddo
             mth_adv(NIN,n_lat_aux_grid-88)= tmp(NIN,130)
             write(*,*)"NNN=",nnn
             write(*,*)"tmp=", tmp(NIN,130)
      end if  

     ttpp=0.0_r4
     work1_g = 0.0_r4

!$OMP PARALLEL DO PRIVATE (iblock,K,J,I)
      do iblock= 1, nblocks_clinic
      do k=1,km
      do j=3,jmt-2
      do i=3,imt-2
          ttpp(i,j,k,iblock)=(axmon(i,j,k,nnn,iblock)+aymon(i,j,k,nnn,iblock)+azmon(i,j,k,nnn,iblock)) & 
                             *work_1(i,j,k,iblock)*dzp(k)*tarea(i,j,iblock)
      enddo
      enddo
      enddo
      enddo
!
      do k=1,km
         call gather_global(buffer_r4_global,ttpp(:,:,k,:), master_task,distrb_clinic) 
         if (mytid == 0 ) then
            do j=1,jmt_global
            do i=1,imt_global
               work1_g(i,j)=work1_g(i,j)+ buffer_r4_global(i,j)
         enddo
         enddo
         end if
      enddo
   
!
    if (mytid == 0 ) then
            DO jj=n_lat_aux_grid-89, 1, -1
            DO j=1,jmt_global
            DO I=1,imt_global
               if (tlatd_g(i,j) < lat_aux_edge(jj) .and. tlatd_g(i,j) >= lat_aux_edge(jj+1) ) then
                  mth_adv(NIN,JJ)=mth_adv(NIN,JJ)-work1_g(I,J)
               end if
            ENDDO
            ENDDO
               mth_adv(NIN,JJ)=mth_adv(NIN,JJ+1)+ mth_adv(NIN,JJ)
            ENDDO
    end if

    END IF

    IF(NIN==1)then

     ttpp=0.0_r4
     work1_g = 0.0_r4
!$OMP PARALLEL DO PRIVATE (iblock,K,J,I)
      do iblock = 1, nblocks_clinic
      do k =1,km
      do j= 1, jmt
      do i= 1, imt
#ifdef ISO
         ttpp(i,j,k,iblock)=(dymon_iso(i,j,k,nnn,iblock)+dxmon_iso(i,j,k,nnn,iblock)) &
                            *dzp(k)*tarea(i,j,iblock)*work_1(i,j,k,iblock)
#else
         ttpp(i,j,k,iblock)=(dymon(i,j,k,nnn,iblock)+dxmon(i,j,k,nnn,iblock)) &
                            *dzp(k)*tarea(i,j,iblock)*work_1(i,j,k,iblock)
#endif
      end do
      end do
      end do
      end do
!
      do k=1,km
         call gather_global(buffer_r4_global,ttpp(:,:,k,:), master_task,distrb_clinic) 
         if ( mytid ==0 ) then
         do j=1,jmt_global
         do i=1,imt_global
            work1_g(i,j)=work1_g(i,j)+ buffer_r4_global(i,j)
         enddo
         enddo
         end if
      enddo
!
    if (mytid == 0 ) then
            DO jj=n_lat_aux_grid-2, 1, -1
            DO j=1,jmt_global
            DO I=1,imt_global
               if (tlatd_g(i,j) < lat_aux_edge(jj) .and. tlatd_g(i,j) >= lat_aux_edge(jj+1) ) then
                  mth_dif(NIN,JJ)=mth_dif(NIN,JJ)-work1_g(I,J)
               end if
            ENDDO
            ENDDO
               mth_dif(NIN,JJ)=mth_dif(NIN,JJ+1)+ mth_dif(NIN,JJ)
            ENDDO
    end if

    END IF
!
!
    IF(NIN==1)then

     ttpp=0.0_r4
     work1_g = 0.0_r4
!$OMP PARALLEL DO PRIVATE (iblock,K,J,I)
      do iblock = 1, nblocks_clinic
      do k =1,km
      do j= 1, jmt
      do i= 1, imt
#ifdef ISO
         ttpp(i,j,k,iblock)=(axmon_iso(i,j,k,nnn,iblock)+aymon_iso(i,j,k,nnn,iblock))  &
                            *dzp(k)*tarea(i,j,iblock)*work_1(i,j,k,iblock)
#endif
      end do
      end do
      end do
      end do
!
      do k=1,km
         call gather_global(buffer_r4_global,ttpp(:,:,k,:), master_task,distrb_clinic) 
         if ( mytid ==0 ) then
         do j=1,jmt_global
         do i=1,imt_global
            work1_g(i,j)=work1_g(i,j)+ buffer_r4_global(i,j)
         enddo
         enddo
         end if
      enddo

!
    if (mytid == 0 ) then
            DO jj=n_lat_aux_grid-2, 1, -1
            DO j=1,jmt_global
            DO I=1,imt_global
               if (tlatd_g(i,j) < lat_aux_edge(jj) .and. tlatd_g(i,j) >= lat_aux_edge(jj+1) ) then
                  mth_adv_iso(NIN,JJ)=mth_adv_iso(NIN,JJ)-work1_g(I,J)
               end if
            ENDDO
            ENDDO
               mth_adv_iso(NIN,JJ)=mth_adv_iso(NIN,JJ+1)+ mth_adv_iso(NIN,JJ)
            ENDDO
    end if
!
     ELSE
    tmp = 0.0_r4
    ttpp = 0.0_r4
    work1_g = 0.0_r4

!$OMP PARALLEL DO PRIVATE (iblock,K,J,I)
      do iblock = 1, nblocks_clinic
      do k =1,km
      do j= 1, jmt
      do i= 1, imt
#ifdef ISO         
         ttpp(i,j,k,iblock)=aymon_iso(i,j,k,nnn,iblock)*dzp(k)*dxu(i,j,iblock)*work_1(i,j,k,iblock)
#else
         ttpp(i,j,k,iblock)=aymon(i,j,k,nnn,iblock)*dzp(k)*dxu(i,j,iblock)*work_1(i,j,k,iblock)
#endif
      enddo
      enddo
      enddo
      enddo
      do k=1,km
         call gather_global(buffer_r4_global,ttpp(:,:,k,:), master_task,distrb_clinic)
         if (mytid == 0 ) then
             do i=1,imt_global
                work1_g(i,130)=work1_g(i,130)+ buffer_r4_global(i,130)
             enddo
         end if
      enddo
     
      if (mytid == 0 ) then
          do i=1,imt_global
            tmp(NIN,130)=tmp(NIN,130)-work1_g(i,130)
          enddo
             mth_adv_iso(NIN,n_lat_aux_grid-88)= tmp(NIN,130)
             write(*,*)"tmp=",tmp(NIN,130)
      end if

     ttpp=0.0_r4
     work1_g = 0.0_r4

!$OMP PARALLEL DO PRIVATE (iblock,K,J,I)
      do iblock= 1, nblocks_clinic
      do k=1,km
      do j=1,jmt
      do i=1,imt
#ifdef ISO
         ttpp(i,j,k,iblock)=(axmon_iso(i,j,k,nnn,iblock)+aymon_iso(i,j,k,nnn,iblock)) &
                             *work_1(i,j,k,iblock)*dzp(k)*tarea(i,j,iblock)
#endif
      enddo
      enddo
      enddo
      enddo

      do k=1,km
         call gather_global(buffer_r4_global,ttpp(:,:,k,:), master_task,distrb_clinic)
         if (mytid == 0 ) then
            do j=1,jmt_global
            do i=1,imt_global
               work1_g(i,j)=work1_g(i,j)+ buffer_r4_global(i,j)
         enddo
         enddo
         end if
      enddo 
    
       if (mytid == 0 ) then
            DO jj=n_lat_aux_grid-89, 1, -1
            DO j=1,jmt_global
            DO I=1,imt_global
               if (tlatd_g(i,j) < lat_aux_edge(jj) .and. tlatd_g(i,j) >= lat_aux_edge(jj+1) ) then
                  mth_adv_iso(NIN,JJ)=mth_adv_iso(NIN,JJ)-work1_g(I,J)
               end if
            ENDDO
            ENDDO
               mth_adv_iso(NIN,JJ)=mth_adv_iso(NIN,JJ+1)+ mth_adv_iso(NIN,JJ)
             ENDDO
    end if
   END IF
      

      IF (NNN==1) then
      do j=1,n_lat_aux_grid - 1
         mth_adv(NIN,j)=mth_adv(NIN,j)*D0*CP*1.0E-15
         mth_dif(NIN,j)=mth_dif(NIN,j)*D0*CP*1.0E-15
#ifdef ISO
         mth_adv_iso(NIN,j)=mth_adv_iso(NIN,j)*D0*CP*1.0E-15
#endif
         mth(NIN,j)=mth_adv(NIN,j)+mth_dif(NIN,j)+mth_adv_iso(NIN,j)
      enddo
      ELSE
      do j=1,n_lat_aux_grid - 1
         mth_adv(NIN,j)=mth_adv(NIN,j)*1.0E-3/34.7
         mth_dif(NIN,j)=mth_dif(NIN,j)*1.0E-3/34.7
#ifdef ISO
         mth_adv_iso(NIN,j)=mth_adv_iso(NIN,j)*1.0E-3/34.7
#endif
         mth(NIN,j)=mth_adv(NIN,j)+mth_dif(NIN,j)+mth_adv_iso(NIN,j)
      enddo
      ENDIF
!
!  end if
      ENDDO !End NIN ZWP 2013-10-18

      deallocate (work1_g,work_1 )
      RETURN
      END SUBROUTINE diag_heat_transport

!========================================
      SUBROUTINE init_diagnostics
!========================================
     integer  :: i_copy,j_dim_sh,grid_error, ntt1,ntt2
     real(r8) :: dlat,southern_edge,  np_minus_northern_edge
     real(r8) :: eps_grid
!
     eps_grid= 1.0D-8
     allocate ( TLATD_G(imt_global,jmt_global) )
     WORK = TLAT/degtorad
     call gather_global (TLATD_G, WORK, master_task,distrb_clinic)


     allocate ( ULATD_G(imt_global,jmt_global) )
     WORK = ULAT/degtorad
     call gather_global (ULATD_G, WORK, master_task,distrb_clinic)

     i_copy = 1
     if ( my_task == master_task ) then
       dlat          = c2 * (TLATD_G(i_copy,jmt_global)-ULATD_G(i_copy,jmt_global))
       southern_edge = ULATD_G(i_copy,jmt_global)
     endif
     call broadcast_scalar (dlat,          master_task)
     call broadcast_scalar (southern_edge, master_task)

   if ( trim(horiz_grid_opt) == 'tripole' ) then

     if ( my_task == master_task )  j_dim_sh = count( ULATD_G(i_copy,:) < c0 )
     call broadcast_scalar (j_dim_sh, master_task)

     grid_error = 0
     if ( my_task == master_task ) then
       do j=jmt_global-j_dim_sh+1, jmt_global
         if(any(abs(TLATD_G(:,j)-TLATD_G(i_copy,j)) > eps_grid))then
           grid_error = -1000
         endif
       enddo
     endif
    call broadcast_scalar (grid_error, master_task)

     np_minus_northern_edge = 90.0_r8 + southern_edge
     if ( np_minus_northern_edge >= p5*dlat ) then
       n_lat_aux_grid = nint( np_minus_northern_edge / dlat )
       dlat = np_minus_northern_edge / dble(n_lat_aux_grid)
       n_lat_aux_grid = 2 * j_dim_sh + n_lat_aux_grid
     else
       n_lat_aux_grid = 2 * j_dim_sh
     endif

   endif

   if ( trim(horiz_grid_opt) == 'lat_lon' ) then
     n_lat_aux_grid = jmt_global
     grid_error = 0
     if ( my_task == master_task ) then
       do j=1,n_lat_aux_grid
         if ( any(TLATD_G(:,j) /= TLATD_G(i_copy,j) ) ) grid_error = -1000
       enddo
     endif
     call broadcast_scalar (grid_error, master_task)
   endif
!-----------------------------------------------------------------------
!
!  allocate arrays
!
!-----------------------------------------------------------------------
   allocate ( lat_aux_edge  (n_lat_aux_grid),  &
              lat_aux_center(n_lat_aux_grid-1) )

   if ( my_task == master_task ) then
     if (trim(horiz_grid_opt) == 'tripole') then

       do j = jmt_global - j_dim_sh +1 , jmt_global
         lat_aux_edge(n_lat_aux_grid+j-jmt_global) = ULATD_G(i_copy, j)
       end do
       jj = n_lat_aux_grid
       do j=n_lat_aux_grid+1-2*j_dim_sh, n_lat_aux_grid - j_dim_sh
         lat_aux_edge(j) = - lat_aux_edge(jj)
         jj = jj - 1
       enddo

       do j = jmt_global - j_dim_sh, jmt_global
         lat_aux_center(n_lat_aux_grid-1+j-jmt_global) = TLATD_G(i_copy, j)
       end do
       jj = n_lat_aux_grid-1
       do j=n_lat_aux_grid-2*j_dim_sh+1, n_lat_aux_grid  - j_dim_sh
         lat_aux_center(j) = - lat_aux_center(jj)
         jj = jj - 1
       enddo

       if ( n_lat_aux_grid == 2 * j_dim_sh ) then
            lat_aux_edge(n_lat_aux_grid+1) = 90.0_r8
       else
         do j=1, n_lat_aux_grid-2*j_dim_sh
           lat_aux_edge(j) = lat_aux_edge(n_lat_aux_grid-2*j_dim_sh+1)+(n_lat_aux_grid-2*j_dim_sh-j+1)*dlat
         enddo
         do j=1, n_lat_aux_grid-2*j_dim_sh
           lat_aux_center(j) = lat_aux_edge(n_lat_aux_grid-2*j_dim_sh+1)+(n_lat_aux_grid-2*j_dim_sh-j+P5)*dlat
         enddo
       endif

     else if (trim(horiz_grid_opt) == 'lat_lon') then
       lat_aux_edge(1)               = 90.0_r8
       lat_aux_edge(2:n_lat_aux_grid+1) = ULATD_G(i_copy,:)
       lat_aux_center = TLATD_G(i_copy,:)
     end if

   endif

   call broadcast_array (lat_aux_edge,   master_task)
   call broadcast_array (lat_aux_center, master_task)
!
!   deallocate ( ULATD_G)
!
   allocate (ng_zonal(2,n_lat_aux_grid-1,km+1))
!
  ng_zonal = 0
  if (mytid == 0 ) then
  DO k= 1, KM
  DO jj=n_lat_aux_grid-2, 1, -1
       ntt1=0
       ntt2=0
       DO j=1,jmt_global
       DO I=1,imt_global
         if (tlatd_g(i,j) < lat_aux_edge(jj) .and. tlatd_g(i,j) >= lat_aux_edge(jj+1) ) then
            ntt1 = ntt1 +  1
            if ( k <= kmt_g(i,j) )then
               ng_zonal(1,JJ,K)= ng_zonal(1,jj,k) +  1
            end if
         end if
!
        if ( basin_g (i,j) == 1 .or. basin_g(i,j) == 2) then
           if (tlatd_g(i,j) < lat_aux_edge(jj) .and. tlatd_g(i,j) >= lat_aux_edge(jj+1) ) then
                  ntt2 = ntt2 +  1
               if ( k <= kmt_g(i,j) )then
                  ng_zonal(2,JJ,K)= ng_zonal(2,jj,k) +  1
               end if
            end if
        end if
       ENDDO
       ENDDO
       if (ntt1==0 ) ng_zonal(1,jj,k) = ng_zonal(1, jj+1,k)
       if (ntt2 < 20 ) ng_zonal(2,jj,k) = ng_zonal(2, jj+1,k)
 ENDDO
 ENDDO
 end if
!
    deallocate(kmt_g, basin_g)
!
  END SUBROUTINE init_diagnostics

#endif

end module diag_mod
