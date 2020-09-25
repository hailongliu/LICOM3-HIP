!  CVS: $Id: forc_mod.F90,v 1.7 2003/08/25 07:47:52 lhl Exp $
module forc_mod
#include <def-undef.h>
!
use precision_mod
use param_mod
!
!
!     ------------------------------------------------------------------
!     Forcing Fields
!     ------------------------------------------------------------------
      real(r8),allocatable,dimension(:,:,:,:) :: su3,sv3,psa3,tsa3,qar3,uva3, &
                   swv3,cld3,sss3,sst3 ,nswv3,dqdt3,chloro3,                  &
                   wspd3,wspdu3,wspdv3,lwv3,seaice3,rain3,snow3,runoff3
!lhl090730
      real(r8),dimension(imt,jmt,max_blocks_clinic)::tsa,uva,qar,cld,&
                    ddd,qqq,sst,dqdt,chloro,lwv,seaice,rain,snow,fresh,runoff, &
                    lthf,sshf !only for output
      real(r8),dimension(imt,jmt,max_blocks_clinic)::BUOYTUR,BUOYSOL
!
      real(r8),allocatable,dimension(:,:,:)::su3_io,sv3_io,psa3_io,tsa3_io,qar3_io,uva3_io, &
                    swv3_io,cld3_io,sss3_io,sst3_io,nswv3_io,dqdt3_io,chloro3_io, &
                    wspdu3_io,wspdv3_io,lwv3_io,seaice3_io,rain3_io,snow3_io,runoff3_io
!
       real(r8),allocatable,dimension(:,:,:,:,:)::restore
       real(r8),allocatable,dimension(:,:,:):: tsf,ssf,su,sv,sss,swv,nswv,USTAR,psa

#if ( defined TIDEMIX )
      real(r8),allocatable,dimension(:,:,:):: wave_dis
#endif
! avoid > 2GB common for 5km version
      real(r8),allocatable,dimension(:,:,:) :: t10,u10,v10,slp,q10,swhf,lwhf
      real(r8),allocatable,dimension(:,:,:) :: precr,precs,rf,si
      real(r8),allocatable,dimension(:,:,:,:) :: buffer3d !ssave
end module forc_mod
