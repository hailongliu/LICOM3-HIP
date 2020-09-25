!-----------------------------------------------------------------------------
!   Processing some variables from flux coupler.
!-----------------------------------------------------------------------------
!
!        By Yongqiang Yu , 16 Apr., 1999
!
!
!
      SUBROUTINE post_cpl

#include <def-undef.h>


use param_mod
use pconst_mod
use tracer_mod
use forc_mod
use buf_mod
use control_licom_mod
use shr_sys_mod
use shr_const_mod
use output_mod,only:spval
use domain
use grid
use blocks
use POP_HaloMod
use POP_GridHorzMod
use distribution
use gather_scatter
use work_mod, only : buffer_real4
use constant_mod, only : boundary_restore
!
      implicit none
      real(r8),dimension(:,:,:),allocatable::tmp_su,tmp_sv
      real(r8) :: ek0
      integer :: iblock, ErrorCode
!
    type (block) :: this_block          ! block information for current block

         allocate(tmp_su(imt,jmt,max_blocks_clinic))
         allocate(tmp_sv(imt,jmt,max_blocks_clinic))
!
        tsf=0.0_r8
        ssf=0.0_r8
        swv=0.0_r8
        tmp_su=0.0_r8
        tmp_sv=0.0_r8

     n=0
!$OMP PARALLEL DO PRIVATE (iblock,this_block,j,i)
    do iblock = 1, nblocks_clinic
        this_block = get_block(blocks_clinic(iblock),iblock)
        do j=this_block%je, this_block%jb, -1
        do i=this_block%ib, this_block%ie
           TSF(i,j,iblock) = (lat1(i,j,iblock)+sen(i,j,iblock)+lwup(i,j,iblock)+ & 
                             lwdn(i,j,iblock)+ netsw(i,j,iblock)+melth(i,j,iblock)-  &
                             iceoff(i,j,iblock)*SHR_CONST_LATICE- & 
                             snow1(i,j,iblock)*SHR_CONST_LATICE) *OD0CP  ! net   heat flux
           SWV(i,j,iblock) =  netsw(i,j,iblock)                          ! net solar radiation
           NSWV(i,j,iblock) = (lat1(i,j,iblock)+sen(i,j,iblock)+lwup(i,j,iblock)+ & 
                             lwdn(i,j,iblock)+ melth(i,j,iblock)-  &
                             iceoff(i,j,iblock)*SHR_CONST_LATICE- & 
                             snow1(i,j,iblock)*SHR_CONST_LATICE) *OD0CP  ! net   heat flux
        if ( boundary_restore == 2) then
           roff(i,j,iblock) = 0.0
           iceoff(i,j,iblock) = 0.0
           SSF(i,j,iblock) =  -(prec(i,j,iblock)+evap(i,j,iblock)+&
                             meltw(i,j,iblock)+roff(i,j,iblock)+     &
                            iceoff(i,j,iblock))*OD0*34.7*1.0D-3/DZP(1)+ &
                            salt(i,j,iblock)*OD0/DZP(1)   ! P+E+melting !linpf 25->DZP(1)
        else
           SSF(i,j,iblock) =  -(prec(i,j,iblock)+evap(i,j,iblock)+&
                             meltw(i,j,iblock)+roff(i,j,iblock)+     &
                            iceoff(i,j,iblock))*OD0*34.7*1.0D-3/DZP(1)+ &
                            salt(i,j,iblock)*OD0/DZP(1)   ! P+E+melting !linpf 25->DZP(1)
        end if
!          SSF(i,j,iblock) =  -(prec(i,j,iblock)+evap(i,j,iblock)+&
!                            meltw(i,j,iblock)+roff(i,j,iblock)+     &
!                           iceoff(i,j,iblock)+salt(i,j,iblock)) &
!                           *34.7*1.0e-3/DZP(1)*OD0   ! P+E+melting !linpf 25->DZP(1)
           tmp_su(i,j,iblock)= taux(i,j,iblock)*cos(anglet(i,j,iblock))+tauy(i,j,iblock)*sin(anglet(i,j,iblock))
           tmp_sv(i,j,iblock)=-tauy(i,j,iblock)*cos(anglet(i,j,iblock))+taux(i,j,iblock)*sin(anglet(i,j,iblock))
        end do
        end do
     end do
        
!!
!!      Update boundary here !!!!
!!
       call POP_HaloUpdate(tsf(:,:,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
       call POP_HaloUpdate(swv(:,:,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
       call POP_HaloUpdate(nswv(:,:,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
       call POP_HaloUpdate(ssf(:,:,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
       call POP_HaloUpdate(tmp_su(:,:,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
       call POP_HaloUpdate(tmp_sv(:,:,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
#ifdef USE_OCN_CARBON
       call POP_HaloUpdate(pco2(:,:,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
#endif         
!!
            lthf = lat1 !latent flux
            sshf = sen !sensible flux
            lwv = lwup + lwdn   !long wave flux
            fresh = ssf
            runoff=roff+iceoff
!

!!linpf 2012Jul26 !2012Jul28
        where(vit(:,:,1,:)<0.5) tsf=spval
        where(vit(:,:,1,:)<0.5) swv=spval
        where(vit(:,:,1,:)<0.5) nswv=spval
        where(vit(:,:,1,:)<0.5) ssf=spval
        where(vit(:,:,1,:)<0.5) lwv=spval
        where(vit(:,:,1,:)<0.5) sshf=spval
        where(vit(:,:,1,:)<0.5) fresh=spval
        where(vit(:,:,1,:)<0.5) runoff=spval
        where(vit(:,:,1,:)<0.5) lthf=spval
!
!
!$OMP PARALLEL DO PRIVATE (iblock)
        do iblock =1 , nblocks_clinic
            call tgrid_to_ugrid(su(:,:,iblock), tmp_su(:,:,iblock),iblock)
            call tgrid_to_ugrid(sv(:,:,iblock), tmp_sv(:,:,iblock),iblock)
        end do
!
!      Update boundary here !!!!
       call POP_HaloUpdate(su(:,:,:) , POP_haloClinic, POP_gridHorzLocSwcorner,&
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
       call POP_HaloUpdate(sv(:,:,:) , POP_haloClinic, POP_gridHorzLocSwcorner,&
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
!
        ! lihuimin, 2012.7.23, ft. lpf
!!linpf 2012Jul29
        where(viv(:,:,1,:)<0.5) su=spval
        where(viv(:,:,1,:)<0.5) sv=spval
!linpf 2012Jul29
!calculate USTAR

!$OMP PARALLEL DO PRIVATE (iblock,j,i)
    do iblock = 1, nblocks_clinic
       DO J = 1, jmt
         DO I = 1,imt
          USTAR(I,J,iblock)=sqrt(sqrt(tmp_su(i,j,iblock)*tmp_su(i,j,iblock)+tmp_sv(i,j,iblock)*tmp_sv(i,j,iblock))*OD0)*vit(i,j,1,iblock) 
         END DO
      END DO        
   end do   
!
!     if (mytid == 5 ) then
!        write(140+mytid,*) ((su(i,j,1),i=3,imt-2),j=3,jmt-2)
!        write(140+mytid,*) "OK"
!        write(140+mytid,*) ((sv(i,j,1),i=3,imt-2),j=3,jmt-2)
!        close(140+mytid)
!     end if
   licomqice = 0.0_r8
   psa       = 0.0_r8
!        if (mytid == 97) then
!          allocate(buffer_real4(imt,jmt))
!          open(35,file='FORCING.DATA', form='unformatted', access='direct', &
!                  recl=(imt-4)*(jmt-4)*4)
!          buffer_real4(:,:) = su(:,:,1)
!          write(35,rec=1) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!          buffer_real4(:,:) = sv(:,:,1)
!          write(35,rec=2) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!          buffer_real4(:,:) = tmp_su(:,:,1)
!          write(35,rec=3) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!          buffer_real4(:,:) = tmp_sv(:,:,1)
!          write(35,rec=4) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!          buffer_real4(:,:) = taux(:,:,1)
!          write(35,rec=5) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!          buffer_real4(:,:) = tauy(:,:,1)
!          write(35,rec=6) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!          buffer_real4(:,:) = lat1(:,:,1)
!          write(35,rec=7) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!          buffer_real4(:,:) = sen(:,:,1)
!          write(35,rec=8) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!          buffer_real4(:,:) = swv(:,:,1)
!          write(35,rec=9) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!          close(35)
!          deallocate(buffer_real4)
!        end if

         deallocate (tmp_su,tmp_sv) !LPF 20120818
!       su = 0.0_r8
!       sv = 0.0_r8
!       tsf= 0.0_r8
!       swv= 0.0_r8
!       nswv= 0.0_r8
!       ssf= 0.0_r8
        return
        end
