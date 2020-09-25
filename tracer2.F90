!  CVS: $Id: tracer.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE TRACER_TDEL(n2)
!     =================
#include <def-undef.h>
use precision_mod 
use param_mod
use pconst_mod
use tracer_mod
use work_mod
use dyn_mod
use isopyc_mod
use forc_mod
use pmix_mod
use msg_mod
use smuvh
use advection
use blocks
use domain
use LICOM_Error_mod
use gather_scatter
use distribution
use buf_mod
#ifdef BIHAR
use hmix_del4
#else
use hmix_del2
#endif
      IMPLICIT NONE
 
      integer     :: n2, iblock,kt, irec
      REAL(r8)    :: AIDIF,C2DTTS,AA,FAW,FIW,ALF,RNCC,ABC,fil_lat1,fil_lat2
      REAL(r8)    :: HDTK(imt,jmt), adv_tt(imt,jmt,km),ek0 ,ek1
      REAL (r8)   :: DT2K(imt,jmt)
#ifdef SSSNORM
!lhl20130913
      REAL(r8)    :: ERR_norm1,ERR_norm2
!lhl20130913
      real (r8),dimension(:,:,:),allocatable::temp11
      real (r8),dimension(:,:),allocatable::temp12 
#endif

!Xiao Chan (Hereinafter XC for short)
      real(r8)    :: LAMDA(imt,jmt,km,max_blocks_clinic),wt1,wt2,adv_y,adv_x,adv_z,adv_x1,adv_x2
      real(r8)    :: LAMDA1(km)
      type (block) :: this_block          ! block information for current block
#ifdef LPFDIAG
      REAL(r8)    :: outsum 
#endif   
      real(r8),dimension(imt,jmt,km,max_blocks_clinic)::VTL_ori !for output dt diffusion
!XC
 
#if (defined BIHAR)
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block,K,hdtk)
   do iblock = 1, nblocks_clinic
   this_block = get_block(blocks_clinic(iblock),iblock)
   do k =1, km
!20200221 LPF-yyq
       call hdifft_del4(k,DT2K,HDTK,ATB(:,:,k,1,iblock),this_block)
       !call hdifft_del4(k,HDTK,ATB(:,:,k,n,iblock),this_block)
!20200221 LPF-yyq
      do j= 3, jmt-2
      do i= 3, imt-2
           TF (I,J,K,IBLOCK) = TF (I,J,K,IBLOCK) + HDTK(I,J)
           dx(i,j,k,N,iblock)=HDTK(i,j) !LPF20160822
      end do
      end do
           !dx(:,:,k,N,iblock)=HDTK(:,:) !LPF20160822
   end do
   end do
!
#else
!
!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block,K,hdtk)
   do iblock = 1, nblocks_clinic
   this_block = get_block(blocks_clinic(iblock),iblock)
   do k =1, km
      call hdifft_del2(k,HDTK,ATB(:,:,k,1,iblock),this_block)
      do j= 3, jmt-2
      do i= 3, imt-2
           TF (I,J,K,IBLOCK) = TF (I,J,K,IBLOCK) + HDTK(I,J)
           dx(i,j,k,N,iblock)=HDTK(i,j) !LPF20160822
      end do
      end do
           !dx(:,:,k,N,iblock)=HDTK(:,:) !LPF20160822
   end do
   end do
!
#endif
      RETURN
      END SUBROUTINE TRACER_TDEL
 
