!  CVS: $Id: readyt.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE setidealvalue
!     =================
!     PREPARATION OF BAROTROPIC AND BAROCLINIC INTEGRATION
 
#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use tracer_mod
use domain
use grid
use blocks
use constant_mod
use operators

!
      IMPLICIT NONE
      REAL(r8)   :: ABCD1
      integer :: iblock
      type (block) :: this_block          ! block information for current block
 
      if (ist == 0 .and. iday ==1) then
        do iblock = 1, nblocks_clinic
        do k=1,km
        do j=1 , jmt
        do i = 1,imt
!           AT(i,j)=28.0*cos(tlat(i,j,1))*
            if(k<10) then
            at(i,j,k,1,iblock) = 0.84*2*to(k)*cos(tlat(i,j,iblock))
            else
             at(i,j,k,1,iblock) = 1.2*to(k)*cos(tlat(i,j,iblock))
            endif
           !at(i,j,k,1,iblock) = to(k)
           at(i,j,k,2,iblock) = 0.2*so(k)
           !at(i,j,k,2,iblock) = so(k)
           !atb(i,j,k,1,iblock) = to(k)
           !atb(i,j,k,1,iblock) = 0.84*2*to(k)*(tlat(i,j,iblock))
            if(k<10) then
            at(i,j,k,1,iblock) = 0.84*2*to(k)*cos(tlat(i,j,iblock))
            else
             at(i,j,k,1,iblock) = 1.2*to(k)*cos(tlat(i,j,iblock))
            endif
           atb(i,j,k,2,iblock) = 0.2*so(k)
           !atb(i,j,k,2,iblock) = so(k)
        end do
        end do
        end do
        end do
      end if

      RETURN
      END SUBROUTINE setidealvalue
 
 
