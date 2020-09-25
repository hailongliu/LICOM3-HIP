!  CVS: $Id: accumm.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE DAILYAVG
!     =================
 
#include <def-undef.h>
use param_mod
use dyn_mod
use tracer_mod
use output_mod
use pconst_mod
use isopyc_mod
#if ( defined TIDEMIX )
use forc_mod, only: su,sv,lthf,sshf,lwv,swv,fresh,runoff,seaice,wave_dis
#else
use forc_mod, only: su,sv,lthf,sshf,lwv,swv,fresh,runoff,seaice
#endif
use domain
use buf_mod,only:prec,evap,iceoff,roff
!use grid, only: FCOR, FCORT 
      IMPLICIT NONE
      integer :: iblock
      real :: Num_total 

    if (daily_accum==.true.) then
        Num_total=float(Num_outputacc)
        if(mytid==0) write(16,*)'num_total=',num_total,Num_outputacc
     !    Num_total=float(num_cpl)
    endif

#if (defined DAILYACC)
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   do iblock = 1, nblocks_clinic
      DO J = 1,JMT
         DO I = 1,IMT

           Z0daily (I,J,iblock)= Z0daily (I,J,iblock)/Num_total
           sudaily (I,J,iblock)= sudaily (I,J,iblock)/Num_total        !U windstress
           svdaily (I,J,iblock)= svdaily (I,J,iblock)/Num_total        !V windstress
           lthfdaily (I,J,iblock)= lthfdaily (I,J,iblock)/Num_total !latent flux
           sshfdaily (I,J,iblock)= sshfdaily (I,J,iblock)/Num_total !sensible flux
           lwvdaily (I,J,iblock)= lwvdaily (I,J,iblock)/Num_total    !long wave flux
           swvdaily (I,J,iblock)= swvdaily (I,J,iblock)/Num_total    !shortwave flux
!       freshdaily(I,J,iblock)= freshdaily (I,J,iblock) + fresh (I,J,iblock)/34.7*1.0e3*DZP(1)/OD0/Num_total !virtual slat flux
!      runoffdaily (I,J,iblock)= runoffdaily (I,J,iblock)/Num_total    !shortwave flux
           mlddaily (I,J,iblock)= mlddaily (I,J,iblock)/Num_total
!           ifracdaily (I,J,iblock)= ifracdaily (I,J,iblock)/Num_total
            precdaily(i,j,iblock)=precdaily(i,j,iblock)/Num_total
            evapdaily(i,j,iblock)=evapdaily(i,j,iblock)/Num_total
            roffdaily(i,j,iblock)=roffdaily(i,j,iblock)/Num_total !+iceoff(i,j,iblock)/Num_total
             
        DO K = 1,KM
              TSdaily (I,J,K,iblock)= TSdaily (I,J,K,iblock)/Num_total
              SSdaily (I,J,K,iblock)= SSdaily (I,J,K,iblock)/Num_total
              USdaily (I,J,K,iblock)= USdaily (I,J,K,iblock)/Num_total
              VSdaily (I,J,K,iblock)= VSdaily (I,J,K,iblock)/Num_total
              WSdaily (I,J,K,iblock)= WSdaily (I,J,K,iblock)/Num_total
        ENDDO

         END DO
      END DO
     END DO
 
#if (defined DAILYBUGDET)

!call POP_HaloUpdate(net , POP_haloClinic, POP_gridHorzLocCenter,&
!                POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,K,N)
   do iblock = 1, nblocks_clinic
        DO J = 1,JMT
          DO I = 1,IMT
           DO K = 1,KM
               pendaily (I,J,K,IBLOCK)= pendaily (I,J,K,IBLOCK)/Num_total !horizontal diffusion
            DO N = 1,NTRA
               tenddaily (I,J,K,N,IBLOCK)= tenddaily (I,J,K,N,IBLOCK)/Num_total !tendency
               axdaily (I,J,K,N,IBLOCK)= axdaily(I,J,K,N,IBLOCK)/Num_total !zonal adv
               aydaily (I,J,K,N,IBLOCK)= aydaily(I,J,K,N,IBLOCK)/Num_total !meridional adv
               azdaily (I,J,K,N,IBLOCK)= azdaily(I,J,K,N,IBLOCK)/Num_total !vertical adv
               dt_convdaily (I,J,K,N,IBLOCK)= dt_convdaily(I,J,K,N,IBLOCK)/Num_total !vertical adv
!LPF20160819
               dxdaily (I,J,K,N,IBLOCK)= dxdaily (I,J,K,N,IBLOCK)/Num_total !horizontal diffusion
               dzdaily (I,J,K,N,IBLOCK)= dzdaily (I,J,K,N,IBLOCK)/Num_total !vertical diffusion
               dt_diffdaily (I,J,K,N,IBLOCK)= dt_diffdaily (I,J,K,N,IBLOCK)/Num_total !vertical diffusion due to dt
           if ( simple_assm ) then
            dt_restoredaily(I,J,K,N,IBLOCK)= dt_restoredaily(I,J,K,N,IBLOCK)/Num_total ! due to assim 
           endif
               if(K==1) net_daily (I,J,N,IBLOCK)= net_daily (I,J,N,IBLOCK)/Num_total !horizontal diffusion
            END DO !N
         END DO !K
       END DO !I
       END DO !J
      END DO !iblock
#endif
!
#endif
      RETURN
      END SUBROUTINE DAILYAVG
 
 
