!  CVS: $Id: tracer_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module tracer_mod
#include <def-undef.h>
!
use precision_mod
use param_mod
!     ------------------------------------------------------------------
!     U V T S H0 W RHO
!     ------------------------------------------------------------------

implicit none

      real(r8),allocatable,dimension(:,:,:,:,:)::atb
      real(r8),allocatable,dimension(:,:,:,:)::net
!mohr
      real(r8),allocatable,dimension(:,:,:,:,:)::at
      real(r8),allocatable,dimension(:,:,:,:,:)::restore_at !LPF20170729
!      real(r8),allocatable,dimension(:,:,:,:,:)::atzwp !ZWP20170218
!
      real(r8),dimension(imt,jmt,km,max_blocks_clinic)::pdensity
      real(r8),allocatable,dimension(:,:,:)::amld
!lhl1204
!
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::tend
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::ax,ay,az
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dx,dy,dz
      real(r8),dimension(imt,jmt,km,max_blocks_clinic)::penetrate
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dt_diff
!
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::ddy
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dt_conv
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dt_restore !20170802
!
#ifdef ISO
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::aay_iso,ddy_iso
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::ax_iso,ay_iso,az_iso
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dx_iso,dy_iso,dz_iso
#endif
!
!     ------------------------------------------------------------------
!     Sea Ice
!     ------------------------------------------------------------------
!
      real(r8),allocatable,dimension(:,:,:),public:: licomqice
      real(r8) FW_norm2
 !
end module tracer_mod
