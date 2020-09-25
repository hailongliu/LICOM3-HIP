!  CVS: $Id: accumm.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE DAILYACCUMM
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
      real :: Num_op2 

    if (daily_accum==.true.) then
        Num_op2= 1 !float(Num_dailyoutput)
     !  if(mytid==0) write(16,*)'num_op2=',num_op2,Num_outputacc
     !    Num_op2=float(num_cpl)
    endif

#if (defined DAILYACC)
!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
   do iblock = 1, nblocks_clinic
      DO J = 1,JMT
         DO I = 1,IMT

           Z0daily (I,J,iblock)= Z0daily (I,J,iblock) + H0(I,J,iblock)/Num_op2
           sudaily (I,J,iblock)= sudaily (I,J,iblock) + su(I,J,iblock)/Num_op2        !U windstress
           svdaily (I,J,iblock)= svdaily (I,J,iblock) + sv(I,J,iblock)/Num_op2        !V windstress
           lthfdaily (I,J,iblock)= lthfdaily (I,J,iblock) + lthf (I,J,iblock)/Num_op2 !latent flux
           sshfdaily (I,J,iblock)= sshfdaily (I,J,iblock) + sshf (I,J,iblock)/Num_op2 !sensible flux
           lwvdaily (I,J,iblock)= lwvdaily (I,J,iblock) + lwv (I,J,iblock)/Num_op2    !long wave flux
           swvdaily (I,J,iblock)= swvdaily (I,J,iblock) + swv (I,J,iblock)/Num_op2    !shortwave flux
!       freshdaily(I,J,iblock)= freshdaily (I,J,iblock) + fresh (I,J,iblock)/34.7*1.0e3*DZP(1)/OD0/Num_op2 !virtual slat flux
!      runoffdaily (I,J,iblock)= runoffdaily (I,J,iblock) + runoff (I,J,iblock)/Num_op2    !shortwave flux
           mlddaily (I,J,iblock)= mlddaily (I,J,iblock) + amld (I,J,iblock)/Num_op2
!           ifracdaily (I,J,iblock)= ifracdaily (I,J,iblock) + seaice(I,J,iblock)/Num_op2
            precdaily(i,j,iblock)=precdaily(i,j,iblock)+prec(i,j,iblock)/Num_op2
            evapdaily(i,j,iblock)=evapdaily(i,j,iblock)+evap(i,j,iblock)/Num_op2
            roffdaily(i,j,iblock)=roffdaily(i,j,iblock)+roff(i,j,iblock)/Num_op2 !+iceoff(i,j,iblock)/Num_op2
             
        DO K = 1,KM
              TSdaily (I,J,K,iblock)= TSdaily (I,J,K,iblock) + AT (I,J,K,1,iblock)*vit(i,j,k,iblock)/Num_op2
              SSdaily (I,J,K,iblock)= SSdaily (I,J,K,iblock) + (AT (I,J,K,2,iblock)*1000.+35.)*vit(i,j,k,iblock)/Num_op2
              USdaily (I,J,K,iblock)= USdaily (I,J,K,iblock) + U (I,J,K,iblock)/Num_op2
              VSdaily (I,J,K,iblock)= VSdaily (I,J,K,iblock) - V(I,J,K,iblock)/Num_op2
              WSdaily (I,J,K,iblock)= WSdaily (I,J,K,iblock) + WS(I,J,K,iblock)/Num_op2
        ENDDO

         END DO
      END DO
     END DO
 
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
!call POP_HaloUpdate(net , POP_haloClinic, POP_gridHorzLocCenter,&
!                POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,K,N)
   do iblock = 1, nblocks_clinic
        DO J = 1,JMT
          DO I = 1,IMT
           DO K = 1,KM
               pendaily (I,J,K,IBLOCK)= pendaily (I,J,K,IBLOCK)+penetrate(i,j,k,IBLOCK)/Num_op2 !horizontal diffusion
            DO N = 1,NTRA
               tenddaily (I,J,K,N,IBLOCK)= tenddaily (I,J,K,N,IBLOCK)+tend(i,j,k,n,IBLOCK)/Num_op2 !tendency
               axdaily (I,J,K,N,IBLOCK)= axdaily(I,J,K,N,IBLOCK)+ax(i,j,k,n,IBLOCK)/Num_op2 !zonal adv
               aydaily (I,J,K,N,IBLOCK)= aydaily(I,J,K,N,IBLOCK)+ay(i,j,k,n,IBLOCK)/Num_op2 !meridional adv
               azdaily (I,J,K,N,IBLOCK)= azdaily(I,J,K,N,IBLOCK)+az(i,j,k,n,IBLOCK)/Num_op2 !vertical adv
               dt_convdaily (I,J,K,N,IBLOCK)= dt_convdaily(I,J,K,N,IBLOCK)+dt_conv(i,j,k,n,IBLOCK)/Num_op2 !vertical adv
!LPF20160819
               dxdaily (I,J,K,N,IBLOCK)= dxdaily (I,J,K,N,IBLOCK)+dx(i,j,k,n,IBLOCK)/Num_op2 !horizontal diffusion
               dzdaily (I,J,K,N,IBLOCK)= dzdaily (I,J,K,N,IBLOCK)+dz(i,j,k,n,IBLOCK)/Num_op2 !vertical diffusion
               dt_diffdaily (I,J,K,N,IBLOCK)= dt_diffdaily (I,J,K,N,IBLOCK)+dt_diff(i,j,k,n,IBLOCK)/Num_op2 !vertical diffusion due to dt
           if ( simple_assm ) then
            dt_restoredaily(I,J,K,N,IBLOCK)= dt_restoredaily(I,J,K,N,IBLOCK)+dt_restore(i,j,k,n,IBLOCK)/Num_op2 ! due to assim 
           endif
               !if(N==1) pendaily (I,J,K,IBLOCK)= pendaily (I,J,K,IBLOCK)+penetrate(i,j,k,IBLOCK)/Num_op2 !horizontal diffusion
               if(K==1) net_daily (I,J,N,IBLOCK)= net_daily (I,J,N,IBLOCK)+net(i,j,N,IBLOCK)/Num_op2 !horizontal diffusion
            END DO !N
         END DO !K
       END DO !I
       END DO !J
      END DO !iblock
#endif
!
#endif
      RETURN
      END SUBROUTINE DAILYACCUMM
 
 
