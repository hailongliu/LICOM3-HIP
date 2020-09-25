!
!-----------------------------------------------------------------------------
!   Preparing oceanic variables for communication between OGCM and coupler.
!-----------------------------------------------------------------------------
!
!        By Yongqiang Yu , 16 Apr., 1999
!

module fluxcpl

#include <def-undef.h>   

      use param_mod
      use pconst_mod
      use buf_mod
      use tracer_mod
      use dyn_mod
      use cdf_mod
      use control_licom_mod
      use msg_mod,only : nproc
      use output_mod,only : qicemon
      use shr_const_mod,only:SHR_CONST_SPVAL !linpf 20120816
      use domain
      use grid
#ifdef USE_OCN_CARBON      
      use coutput_mod, only : uptake
#endif      
!
!
contains

  SUBROUTINE flux_cpl

!
        implicit none
!
    integer :: iblock
    real(r8), allocatable :: tmp_u(:,:), tmp_v(:,:)
!
    allocate (tmp_u(imt,jmt),tmp_v(imt,jmt))
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
    do iblock = 1, nblocks_clinic
        do j=1,jmt
        do i=1,imt
           if (licomqice(i,j,iblock) .gt. 0.0) then
               q(i,j,iblock)= licomqice(i,j,iblock)*D0*CP*DZP(1)/86400.*vit(i,j,1,iblock)*dble(num_cpl)
           else
               q(i,j,iblock)= (tbice-at(i,j,1,1,iblock))*D0*CP*DZP(1)/86400.*vit(i,j,1,iblock)*dble(num_cpl)
           endif
        end do
        end do
    end do
!
        T_CPL = (AT(:,:,1,1,:) + 273.15)*vit(:,:,1,:)
        S_CPL = (AT(:,:,1,2,:)*1000. + 35.)*vit(:,:,1,:)
        U_CPL = 0.
        V_CPL = 0.
!linpf 2012Jul26
 
!$OMP PARALLEL DO PRIVATE (IBLOCK,tmp_u,tmp_v)
    do iblock = 1, nblocks_clinic
      call ugrid_to_tgrid(tmp_u,u(:,:,1,iblock),iblock,1)
      call ugrid_to_tgrid(tmp_v,v(:,:,1,iblock),iblock,1)
      do j = 1, jmt
      do i = 1, imt
         u_cpl(i,j,iblock) = tmp_u(i,j)*cos(anglet(i,j,iblock)) + tmp_v(i,j)*sin(anglet(i,j,iblock))
         v_cpl(i,j,iblock) = tmp_u(i,j)*sin(anglet(i,j,iblock)) - tmp_v(i,j)*cos(anglet(i,j,iblock))
      end do
      end do
    end do
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,tmp_u,tmp_v)
    do iblock = 1, nblocks_clinic
        do j= 2, jmt-1
        do i= 2, imt-1
           tmp_u(i,j)  =   (h0(i+1,j,iblock)-h0(i-1,j,iblock)) /(hun(i,j,iblock)+hun(i+1,j,iblock))
           tmp_v(i,j)  =   (h0(i,j+1,iblock)-h0(i,j-1,iblock)) /(hue(i,j,iblock)+hue(i,j-1,iblock))
           dhdx(i,j,iblock) = tmp_u(i,j)*cos(anglet(i,j,iblock)) + tmp_v(i,j)*sin(anglet(i,j,iblock))
           dhdy(i,j,iblock) = tmp_u(i,j)*sin(anglet(i,j,iblock)) - tmp_v(i,j)*cos(anglet(i,j,iblock))
        end do
        end do
     end do

        ! TODO, consider here
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
    do iblock = 1,nblocks_clinic
        do j=1,jmt
         do i=1,imt
           if(vit(i,j,1,iblock)<0.5) dhdx(i,j,iblock)=0.0_r8 !SHR_CONST_SPVAL  
           if(vit(i,j,1,iblock)<0.5) dhdy(i,j,iblock)=0.0_r8 !SHR_CONST_SPVAL
         enddo
        enddo
    enddo
!
!
     where (q(:,:,:) > 0 ) qicemon(:,:,:) = qicemon(:,:,:)+ q(:,:,:)
     licomqice = 0.0
     deallocate(tmp_u,tmp_v)
        
#ifdef USE_OCN_CARBON
      co2_cpl(:,:) = uptake(:,:)
#endif          
        return

  END subroutine flux_cpl


end module fluxcpl

