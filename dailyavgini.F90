!  CVS: $Id: mm00.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ===================
      SUBROUTINE DAILYAVGINI 
!     ===================
 
#include <def-undef.h>
use param_mod
use pconst_mod
use constant_mod
use output_mod
      IMPLICIT NONE

#if (defined DAILYACC)

               z0daily = 0.0
               mlddaily = 0.0
               sudaily = 0.0
               svdaily = 0.0
               lthfdaily = 0.0
               sshfdaily = 0.0
               lwvdaily = 0.0
               swvdaily = 0.0
               TSdaily = 0.0
               SSdaily = 0.0
               USdaily = 0.0
               VSdaily = 0.0
               WSdaily = 0.0
               precdaily = 0.0
               evapdaily = 0.0
               roffdaily = 0.0
               
#if (defined DAILYBUGDET)

     if (.not. allocated(tenddaily)) then 
      allocate(tenddaily(imt,jmt,km,ntra,max_blocks_clinic))
      tenddaily=0.0
     endif
     if (.not. allocated(dt_diffdaily)) allocate(dt_diffdaily(imt,jmt,km,ntra,max_blocks_clinic))
     if (.not. allocated(axdaily)) allocate(axdaily(imt,jmt,km,ntra,max_blocks_clinic))
     if (.not. allocated(aydaily)) allocate(aydaily(imt,jmt,km,ntra,max_blocks_clinic))
     if (.not. allocated(azdaily)) allocate(azdaily(imt,jmt,km,ntra,max_blocks_clinic))
     if (.not. allocated(dxdaily)) allocate(dxdaily(imt,jmt,km,ntra,max_blocks_clinic))
     if (.not. allocated(dzdaily)) allocate(dzdaily(imt,jmt,km,ntra,max_blocks_clinic))
     if (.not. allocated(dt_convdaily)) allocate(dt_convdaily(imt,jmt,km,ntra,max_blocks_clinic))
     if (.not. allocated(pendaily)) allocate(pendaily(imt,jmt,km,max_blocks_clinic))
     if (.not. allocated(net_daily)) then
       allocate(net_daily(imt,jmt,ntra,max_blocks_clinic))
      net_daily=0.0
     endif
     if (.not. allocated(dt_restoredaily)) then 
      allocate(dt_restoredaily(imt,jmt,km,ntra,max_blocks_clinic))
      dt_restoredaily=0.0
     endif   

               tenddaily=0.0
               dt_diffdaily=0.0
               axdaily=0.0
               aydaily=0.0
               azdaily=0.0
               dxdaily=0.0
               dzdaily=0.0
               pendaily=0.0
               dt_convdaily=0.0
               if ( simple_assm ) dt_restoredaily=0.0
               !if ( simple_assm==.true.) dt_restoredaily=0.0
               net_daily=0.0 
#endif

#endif

      RETURN
      END SUBROUTINE DAILYAVGINI
 
 
