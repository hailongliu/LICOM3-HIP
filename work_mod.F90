!  CVS: $Id: work_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module work_mod
#include <def-undef.h>
!
use precision_mod
use param_mod
!
!
!     ------------------------------------------------------------------
!     Working Arrays
!     ------------------------------------------------------------------
      real(r8),dimension(imt,jmt,max_blocks_clinic):: PXB,PYB,PAX,PAY,WHX,WHY,WGP
      real(r8),dimension(:,:,:,:),allocatable::wka
      real(r8),dimension(:,:,:,:),allocatable:: work_1,work_2,work_3,temp
      real(r8),dimension(:,:,:,:),allocatable:: tmp1,tmp2,uk,vk
      real(r8),dimension(:,:,:),allocatable :: work
      real(r8):: wkk(kmp1)
      real(r8),dimension(:,:,:,:),allocatable:: WKB,WKC,WKD,TF
      real(r8),dimension(:,:,:),allocatable::stf
      real(r4),dimension(:,:),allocatable:: buffer_real4
      real(r8),dimension(:,:), allocatable :: work1_g, work2_g, work3_g
!
end module work_mod
