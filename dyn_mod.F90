!  CVS: $Id: dyn_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module dyn_mod
#include <def-undef.h>
use precision_mod
use param_mod
!     ------------------------------------------------------------------
!     U V T S H0 W RHO
!     ------------------------------------------------------------------
      real(r8),allocatable,dimension(:,:,:)::ub,vb,ubp,vbp,h0p
      real(r8),allocatable,dimension(:,:,:,:)::up,vp
      real(r8),allocatable,dimension(:,:,:,:)::ws
      real(r8),allocatable,dimension(:,:,:)::h0l,h0f,h0bl,h0bf,SBCX,BBCX,SBCY,BBCY
      real(r8),allocatable,dimension(:,:,:,:)::utl,vtl,utf,vtf
!


!
      real(r8),allocatable,dimension(:,:) :: buffer
!
      real(r8),allocatable,dimension(:,:,:)::h0
      real(r8),allocatable,dimension(:,:,:,:)::u,v
!
!
!     ------------------------------------------------------------------
!     Pressure gradient
!     ------------------------------------------------------------------
      real(r8),dimension(:,:,:,:),allocatable::gg,dlu,dlv
      real(r8),dimension(:,:,:),allocatable::dlub,dlvb
end module dyn_mod

