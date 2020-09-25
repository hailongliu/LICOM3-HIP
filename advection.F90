!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module advection

!BOP
! !MODULE: advection
!
! !DESCRIPTION:
!  This module contains arrays and variables necessary for performing
!  advection of momentum and tracer quantities.  Currently, the
!  module supports leapfrog centered advection of momentum and
!  both leapfrog centered advection and third-order upwinding of
!  tracers.
!
! !REVISION HISTORY:
!  SVN:$Id: advection.F90 28439 2011-05-18 21:40:58Z njn01 $

! !USES:
   use precision_mod
   use param_mod
   use pconst_mod
   use constant_mod
   use grid
   use blocks
   use tracer_mod, only : ax,ay,az
   use dyn_mod, only : h0
   use LICOM_Error_mod

   implicit none
   private
   save

   public :: advection_momentum, advection_tracer
!

!EOP
!BOC

    contains

!-----------------------------------------------------------------------
!EOC
      subroutine advection_momentum(uuu,vvv,www,adv_uu,adv_vv,iblock)
!
      real(r8) ,intent(in) :: uuu(imt,jmt,km), vvv(imt,jmt,km),www(imt,jmt,km)
      real(r8), intent(out):: adv_uu(imt,jmt,km), adv_vv(imt,jmt,km)
      real(r8) :: u_wface(imt,jmt,km), v_sface(imt,jmt,km)
      real(r8) :: adv_z1, adv_z2,adv_z3,adv_z4
      integer, intent(in)  :: iblock
!
     adv_uu = c0
     adv_vv = c0
!
     if ( adv_momentum(1:8) == 'centered' ) then
      do k=1, km
      do j= 1,jmt-1
      do i= 2,imt
            u_wface(i,j,k) = (uuu(i-1,j,k) + uuu(i,j,k))*P25*hue(i-1,j,iblock)
            v_sface(i,j,k) = (vvv(i,j,k) + vvv(i,j+1,k))*P25*hun(i,j+1,iblock)
      end do
      end do
      end do
     else if ( adv_momentum(1:4) == 'flux' ) then
      do k=1, km
      do j= 1,jmt-1
      do i= 2,imt
            u_wface(i,j,k) = (uuu(i-1,j,k)*dyu(i-1,j,iblock) + uuu(i,j,k)*dyu(i,j,iblock))*P25
            v_sface(i,j,k) = (vvv(i,j,k)*dxu(i,j,iblock) + vvv(i,j+1,k)*dxu(i,j+1,iblock))*P25
      end do
      end do
      end do
     else
        call exit_licom(sigAbort,'The false advection option for tracer')
     end if
!
      if ( adv_momentum(1:8) == 'centered' ) then
        do k = 1,km
        do j= 3,jmt-2
        do i= 3,imt-2
           adv_uu(i,j,k) = (-u_wface(i  ,j,k)*(uuu(i  ,j,k)-uuu(i-1,j,k))                   &
                            -u_wface(i+1,j,k)*(uuu(i+1,j,k)-uuu(i  ,j,k))                   &
                            -v_sface(i,j  ,k)*(uuu(i,j+1,k)-uuu(i,j  ,k))                   &
                            -v_sface(i,j-1,k)*(uuu(i,j  ,k)-uuu(i,j-1,k)))*uarea_r(i,j,iblock)
           adv_vv(i,j,k) = (-u_wface(i  ,j,k)*(vvv(i  ,j,k)-vvv(i-1,j,k))                   &
                            -u_wface(i+1,j,k)*(vvv(i+1,j,k)-vvv(i  ,j,k))                   &
                            -v_sface(i,j  ,k)*(vvv(i,j+1,k)-vvv(i,j  ,k))                   &
                            -v_sface(i,j-1,k)*(vvv(i,j  ,k)-vvv(i,j-1,k)))*uarea_r(i,j,iblock)
!
           if (k==1 )then
               adv_z1=0.0D0
               adv_z3=0.0D0
           else
               adv_z1=www (I,J,K)* (uuu(I,J,K -1) - uuu(I,J,K))
               adv_z3=www (I,J,K)* (vvv(I,J,K -1) - vvv(I,J,K))
           end if
!
           if (k==km )then
                adv_z2=0.0D0
                adv_z4=0.0D0
           else
                adv_z2=www(I,J,K+1)*(uuu(I,J,K ) - uuu(I,J,K+1))
                adv_z4=www(I,J,K+1)*(vvv(I,J,K ) - vvv(I,J,K+1))
           end if
           adv_uu(i,j,k) = adv_uu(i,j,k) - P5*ODZP(K)* (adv_z1+adv_z2)
           adv_vv(i,j,k) = adv_vv(i,j,k) - P5*ODZP(K)* (adv_z3+adv_z4)
        end do
        end do
        end do
      else if ( adv_momentum(1:4) == 'flux' ) then
        do k = 1,km
        do j= 3,jmt-2
        do i= 3,imt-2
           adv_uu(i,j,k) = (-u_wface(i  ,j,k)*(uuu(i  ,j,k)+uuu(i-1,j,k))    &
                            +u_wface(i+1,j,k)*(uuu(i+1,j,k)+uuu(i  ,j,k))    &
                            -v_sface(i,j-1,k)*(uuu(i,  j,k)+uuu(i,j-1,k))    &
                            +v_sface(i,j  ,k)*(uuu(i,j+1,k)+uuu(i,  j,k)))*uarea_r(i,j,iblock)
 
           adv_vv(i,j,k) = (-u_wface(i  ,j,k)*(vvv(i,j  ,k)+vvv(i-1,j,k))    &
                            +u_wface(i+1,j,k)*(vvv(i+1,j,k)+vvv(i  ,j,k))    &
                            -v_sface(i,j-1,k)*(vvv(i,  j,k)+vvv(i,j-1,k))    &
                            +v_sface(i,j  ,k)*(vvv(i,j+1,k)+vvv(i,  j,k)))*uarea_r(i,j,iblock)
 
           if (k ==1 ) then
               adv_z1=0.0D0
               adv_z3=0.0D0
           else
               adv_z1=www (I,J,K)* (uuu(I,J,K -1) + uuu(I,J,K))*P5
               adv_z3=www (I,J,K)* (vvv(I,J,K -1) + vvv(I,J,K))*P5
           end if
!
           if (k == km ) then
                adv_z2=0.0D0
                adv_z4=0.0D0
           else
                adv_z2=www(I,J,K+1)*(uuu(I,J,K ) + uuu(I,J,K+1))*P5
                adv_z4=www(I,J,K+1)*(vvv(I,J,K ) + vvv(I,J,K+1))*P5
           end if
           adv_uu(i,j,k) = adv_uu(i,j,k) - ODZP(K)* (adv_z2-adv_z1)
           adv_vv(i,j,k) = adv_vv(i,j,k) - ODZP(K)* (adv_z4-adv_z3)
        end do
        end do
        end do
      else
        write(6,*) "adv_momentum =", adv_momentum
        call exit_licom(sigAbort,'The false advection option for momentum')
      end if
!
      end subroutine advection_momentum
!
!
!
      !subroutine advection_tracer(uuu,vvv,www,ttt,adv_tt,iblock,mtracer, this_block,ax,ay,az)
      subroutine advection_tracer(uuu,vvv,www,ttt,adv_tt,iblock,mtracer, this_block) !LPF20160811
!
      real(r8) ,intent(in) :: uuu(imt,jmt,km), vvv(imt,jmt,km),www(imt,jmt,km+1),ttt(imt,jmt,km)
      integer, intent(in)  :: iblock,mtracer
      real(r8), intent(out):: adv_tt(imt,jmt,km)
      real(r8) :: u_wface(imt,jmt,km),v_sface(imt,jmt,km)
      real(r8) :: adv_z1, adv_z2, sumv
      real(r8)    :: LAMDA,wt1,wt2,adv_z, ek1,ek0,ek2,sflux(imt,jmt),nflux(imt,jmt)
      real(r8),dimension(:,:,:),allocatable :: adv_xy1,adv_xy2,adv_xy3,adv_xy4, at00
      real(r8),dimension(:,:,:),allocatable :: adv_zz,atz, adv_za,adv_zb1,adv_zb2,atmaxz,atminz
      real(r8),dimension(:,:,:),allocatable :: adv_x0,adv_y0,atmax,atmin,adv_xx,adv_yy
      real(r8),dimension(:,:,:),allocatable :: adv_c1,adv_c2,adv_zc
      type(block) :: this_block
!
     adv_tt=0.0_r8
!
      do k=1, km
      do j=1, jmt
      do i=2, imt-1
         if ( adv_tracer(1:8) == 'centered' .or. adv_tracer(1:5) == 'tspas') then
            v_sface(i,j,k) = (vvv(i,j,k) + vvv(i+1,j,k))*hts(i,j,iblock)*P25
         else if ( adv_tracer(1:4) == 'flux' ) then
            v_sface(i,j,k) = (vvv(i,j,k)*dxu(i,j,iblock) + vvv(i+1,j,k)*dxu(i+1,j,iblock))*P25
         end if
      end do
      end do
      do j= 2,jmt-1
      do i= 1,imt
         if ( adv_tracer(1:8) == 'centered' .or. adv_tracer(1:5) == 'tspas') then
            u_wface(i,j,k) = (uuu(i,j-1,k) + uuu(i,j,k))*htw(i,j,iblock)*P25
         else if ( adv_tracer(1:4) == 'flux' ) then
            u_wface(i,j,k) = (uuu(i,j-1,k)*dyu(i,j-1,iblock) + uuu(i,j,k)*dyu(i,j,iblock))*P25
         end if
      end do
      end do
      end do
!
      if ( adv_tracer(1:8) == 'centered' ) then
        do k = 1,km
        do j= 3,jmt-2
        do i= 3,imt-2
           adv_tt(i,j,k) = (-u_wface(i  ,j,k)*(ttt(i  ,j,k)-ttt(i-1,j,k))                      &
                            -u_wface(i+1,j,k)*(ttt(i+1,j,k)-ttt(i  ,j,k))                      &
                            -v_sface(i,j  ,k)*(ttt(i,j+1,k)-ttt(i,j  ,k))                      &
                            -v_sface(i,j-1,k)*(ttt(i,j  ,k)-ttt(i,j-1,k)))*tarea_r(i,j,iblock)
           ax(i,j,k,mtracer,iblock) =  ax(i,j,k,mtracer,iblock) +  &
                               (-u_wface(i  ,j,k)*(ttt(i  ,j,k)-ttt(i-1,j,k))                      &
                                -u_wface(i+1,j,k)*(ttt(i+1,j,k)-ttt(i  ,j,k)))*tarea_r(i,j,iblock)/dble(NSS)
           ay(i,j,k,mtracer,iblock) =  ay(i,j,k,mtracer,iblock) +  &
                               (-v_sface(i,j  ,k)*(ttt(i,j+1,k)-ttt(i,j  ,k))                      &
                                -v_sface(i,j-1,k)*(ttt(i,j  ,k)-ttt(i,j-1,k)))*tarea_r(i,j,iblock)/dble(NSS)
!
           if (k==1 )then
               adv_z1=0.0D0
           else
               adv_z1=www (I,J,K)* (ttt(I,J,K -1) - ttt(I,J,K))
           end if
!
           if (k==km )then
                adv_z2=0.0D0
           else
                adv_z2=www(I,J,K+1)*(ttt(I,J,K ) - ttt(I,J,K+1))
           end if
!
           adv_tt(i,j,k) = adv_tt(i,j,k) - P5*ODZP(K)* (adv_z1+adv_z2)
           az(i,j,k,mtracer,iblock) =  az(i,j,k,mtracer,iblock) - P5*ODZP(K)* (adv_z1+adv_z2)/dble(NSS)
        end do
        end do
        end do
      else if ( adv_tracer(1:4) == 'flux' ) then
        do k = 1,km
        do j= 3,jmt-2
        do i= 3,imt-2
           adv_tt(i,j,k) = (-u_wface(i  ,j,k)*(ttt(i  ,j,k)+ttt(i-1,j,k))   &
                            +u_wface(i+1,j,k)*(ttt(i+1,j,k)+ttt(i  ,j,k))   &
                            -v_sface(i,j-1,k)*(ttt(i,  j,k)+ttt(i,j-1,k))   &
                            +v_sface(i,j  ,k)*(ttt(i,j+1,k)+ttt(i,  j,k)))*tarea_r(i,j,iblock)
!
           if (k==1 )then
               adv_z1=0.0D0
           else
               adv_z1=www (I,J,K)* (ttt(I,J,K -1) + ttt(I,J,K))*P5
           end if
!
           if (k==km )then
                adv_z2=0.0D0
           else
                adv_z2=www(I,J,K+1)*(ttt(I,J,K ) + ttt(I,J,K+1))*P5
           end if
!
           adv_tt(i,j,k) = adv_tt(i,j,k) - ODZP(K)* (adv_z2-adv_z1)
!
        end do
        end do
        end do
      else if (adv_tracer(1:5) == 'tspas' ) then
!
      allocate ( adv_c1(imt,jmt,km),adv_c2(imt,jmt,km),adv_zc(imt,jmt,km))
      allocate ( adv_x0(imt,jmt,km),adv_y0(imt,jmt,km))
      allocate ( adv_xx(imt,jmt,km), adv_yy(imt,jmt,km))
      allocate ( adv_xy1(imt,jmt,km), adv_xy2(imt,jmt,km), adv_xy3(imt,jmt,km), adv_xy4(imt,jmt,km))
      allocate ( at00(imt,jmt,km), atmax(imt,jmt,km), atmin(imt,jmt,km) )
      allocate (adv_zz(imt,jmt,km), adv_za(imt,jmt,km),adv_zb1(imt,jmt,km) )
      allocate (adv_zb2(imt,jmt,km),atmaxz(imt,jmt,km), atminz(imt,jmt,km), atz(imt,jmt,km))
!
      nflux=0.0D0
      sflux=0.0D0
      DO K = 1,km
        DO J = 2, jmt-1
          DO I = 2, imt-1
            adv_x0(i,j,k)=(ttt(i+1,j,k)+ttt(i,j,k))*u_wface(i+1,j,k)*tarea_r(i,j,iblock) &
                       -(ttt(i,j,k)+ttt(i-1,j,k))*u_wface(i,j,k)*tarea_r(i,j,iblock)
            adv_y0(i,j,k)=((ttt(i,j+1,k)+ttt(i,j,k))*v_sface(i,j,k) &
                       -(ttt(i,j,k)+ttt(i,j-1,k))*v_sface(i,j-1,k))*tarea_r(i,j,iblock)
            adv_xy1(i,j,k)=-dts*(ttt(i+1,j,k)-ttt(i,j,k))*2.0_r8*tarea_r(i,j,iblock)* &
                                  u_wface(i+1,j,k)*u_wface(i+1,j,k)/(htw(i+1,j,iblock)*hun(i+1,j,iblock))
            adv_xy2(i,j,k)= dts*(ttt(i,j,k)-ttt(i-1,j,k))*2.0_r8*tarea_r(i,j,iblock)*  &
                                  u_wface(i,j,k)*u_wface(i,j,k)/(htw(i,j,iblock)*hun(i,j,iblock))
            adv_xy3(i,j,k)=-dts*(ttt(i,j+1,k)-ttt(i,j,k))*2.0_r8*tarea_r(i,j,iblock)* &
                                  v_sface(i,j,k)*v_sface(i,j,k)/(hts(i,j,iblock)*hue(i,j,iblock))
            adv_xy4(i,j,k)= dts*(ttt(i,j,k)-ttt(i,j-1,k))*2.0_r8*tarea_r(i,j,iblock)* &
                                  v_sface(i,j-1,k)*v_sface(i,j-1,k)/(hts(i,j-1,iblock)*hue(i,j-1,iblock))
            adv_c1(i,j,k)=-ttt(i,j,k)*(u_wface(i+1,j,k)-u_wface(i,j,k))*tarea_r(i,j,iblock)*2.0_r8
            adv_c2(i,j,k)=-ttt(i,j,k)*(v_sface(i,j,k)-v_sface(i,j-1,k))*tarea_r(i,j,iblock)*2.0_r8


            if (k==1) then
              adv_za (i,j,k)=-0.5*ODZP(1)*WWW(I,J,2)*(ttt(I,J,2)+ttt(I,J,1))
              adv_zb1(i,j,k)=0
              adv_zb2(i,j,k)= 0.5*ODZP(1)*WWW(I,J,2)*WWW(I,J,2)*ODZT(2) &
                                *(ttt(I,J,1)-ttt(I,J,2))*dts
              adv_zc (i,j,k)=     ODZP(1)*ttt(I,J,1)*WWW(I,J,2)
            else if (k==km) then
              adv_za (i,j,k)= 0.5*ODZP(km)*WWW(I,J,km)*(ttt(I,J,km)+ttt(I,J,km-1))
              adv_zb1(i,j,k)=-0.5*ODZP(km)*WWW(I,J,km)*WWW(I,J,km)*ODZT(km  ) &
                                *(ttt(I,J,km-1)-ttt(I,J,km))*dts
              adv_zb2(i,j,k)=0
              adv_zc (i,j,k)=    -ODZP(km)*ttt(I,J,km)*WWW(I,J,km)
            else
              adv_za (i,j,k)= 0.5*ODZP(k)*WWW(I,J,k  )*(ttt(I,J,k)+ttt(I,J,k-1)) &
                         -0.5*ODZP(k)*WWW(I,J,k+1)*(ttt(I,J,k)+ttt(I,J,k+1))
              adv_zb1(i,j,k)=-0.5*ODZP(k)*WWW(I,J,k  )*WWW(I,J,k  )*ODZT(k  ) &
                                *(ttt(I,J,k-1)-ttt(I,J,k))*dts
              adv_zb2(i,j,k)= 0.5*ODZP(k)*WWW(I,J,k+1)*WWW(I,J,k+1)*ODZT(k+1) &
                                *(ttt(I,J,k)-ttt(I,J,k+1))*dts
              adv_zc (i,j,k)=    -ODZP(k)*ttt(I,J,k)*(WWW(I,J,k)-WWW(I,J,k+1))
            endif
!
            adv_xx(i,j,k)=-(adv_x0(i,j,k)+adv_xy1(i,j,k)+adv_xy2(i,j,k)+adv_c1(i,j,k))
            adv_yy(i,j,k)=-(adv_y0(i,j,k)+adv_xy3(i,j,k)+adv_xy4(i,j,k)+adv_c2(i,j,k))
            adv_zz(i,j,k)=-(adv_zb1(i,j,k)+adv_zb2(i,j,k)+adv_za(i,j,k)+adv_zc(i,j,k))
!           adv_xx(i,j,k)=-(adv_x0(i,j,k)+adv_xy1(i,j,k)+adv_xy2(i,j,k))
!           adv_yy(i,j,k)=-(adv_y0(i,j,k)+adv_xy3(i,j,k)+adv_xy4(i,j,k))
!           adv_zz(i,j,k)=-(adv_zb1(i,j,k)+adv_zb2(i,j,k)+adv_za(i,j,k))
            at00(i,j,k)=ttt(i,j,k)+(adv_xx(i,j,k)+adv_yy(i,j,k)+adv_zz(i,j,k))*dts
          ENDDO
        ENDDO
      ENDDO


      !az(:,:,:,mtracer,iblock)=  adv_xy4 !LPF20160812 wrong just for test
!
      WT1 = -1.0D10
      WT2 = +1.0D10
      DO K = 1,km
        DO J = 2, jmt-1
          DO I = 2,imt-1
           if ( k == 1 ) then
             atmax(i,j,k)=max(ttt(i,j,k)*vit(i,j,k,iblock)+(1.0D0-vit(i,j,k,iblock))*WT1, & 
                              ttt(i,j-1,k)*vit(i,j-1,k,iblock)+(1.0D0-vit(i,j-1,k,iblock))*WT1, & 
                              ttt(i,j+1,k)*vit(i,j+1,k,iblock)+(1.0D0-vit(i,j+1,k,iblock))*WT1, & 
                              ttt(i-1,j,k)*vit(i-1,j,k,iblock)+(1.0D0-vit(i-1,j,k,iblock))*WT1, &
                              ttt(i+1,j,k)*vit(i+1,j,k,iblock)+(1.0D0-vit(i+1,j,k,iblock))*WT1, &
                              ttt(i,j,k+1)*vit(i,j,k+1,iblock)+(1.0D0-vit(i,j,k+1,iblock))*WT1)
             atmin(i,j,k)=min(ttt(i,j,k)*vit(i,j,k,iblock)+(1.0D0-vit(i,j,k,iblock))*WT2, & 
                              ttt(i,j-1,k)*vit(i,j-1,k,iblock)+(1.0D0-vit(i,j-1,k,iblock))*WT2, & 
                              ttt(i,j+1,k)*vit(i,j+1,k,iblock)+(1.0D0-vit(i,j+1,k,iblock))*WT2, & 
                              ttt(i-1,j,k)*vit(i-1,j,k,iblock)+(1.0D0-vit(i-1,j,k,iblock))*WT2, &
                              ttt(i+1,j,k)*vit(i+1,j,k,iblock)+(1.0D0-vit(i+1,j,k,iblock))*WT2, &
                              ttt(i,j,k+1)*vit(i,j,k+1,iblock)+(1.0D0-vit(i,j,k+1,iblock))*WT2)
           else if ( k == km) then
             atmax(i,j,k)=max(ttt(i,j,k)*vit(i,j,k,iblock)+(1.0D0-vit(i,j,k,iblock))*WT1, & 
                              ttt(i,j-1,k)*vit(i,j-1,k,iblock)+(1.0D0-vit(i,j-1,k,iblock))*WT1, & 
                              ttt(i,j+1,k)*vit(i,j+1,k,iblock)+(1.0D0-vit(i,j+1,k,iblock))*WT1, & 
                              ttt(i-1,j,k)*vit(i-1,j,k,iblock)+(1.0D0-vit(i-1,j,k,iblock))*WT1, &
                              ttt(i+1,j,k)*vit(i+1,j,k,iblock)+(1.0D0-vit(i+1,j,k,iblock))*WT1, &
                              ttt(i,j,k-1)*vit(i,j,k-1,iblock)+(1.0D0-vit(i,j,k-1,iblock))*WT1)
             atmin(i,j,k)=min(ttt(i,j,k)*vit(i,j,k,iblock)+(1.0D0-vit(i,j,k,iblock))*WT2, & 
                              ttt(i,j-1,k)*vit(i,j-1,k,iblock)+(1.0D0-vit(i,j-1,k,iblock))*WT2, & 
                              ttt(i,j+1,k)*vit(i,j+1,k,iblock)+(1.0D0-vit(i,j+1,k,iblock))*WT2, & 
                              ttt(i-1,j,k)*vit(i-1,j,k,iblock)+(1.0D0-vit(i-1,j,k,iblock))*WT2, &
                              ttt(i+1,j,k)*vit(i+1,j,k,iblock)+(1.0D0-vit(i+1,j,k,iblock))*WT2, &
                              ttt(i,j,k-1)*vit(i,j,k-1,iblock)+(1.0D0-vit(i,j,k-1,iblock))*WT2 )
           else 
             atmax(i,j,k)=max(ttt(i,j,k)*vit(i,j,k,iblock)+(1.0D0-vit(i,j,k,iblock))*WT1, & 
                              ttt(i,j-1,k)*vit(i,j-1,k,iblock)+(1.0D0-vit(i,j-1,k,iblock))*WT1, & 
                              ttt(i,j+1,k)*vit(i,j+1,k,iblock)+(1.0D0-vit(i,j+1,k,iblock))*WT1, & 
                              ttt(i-1,j,k)*vit(i-1,j,k,iblock)+(1.0D0-vit(i-1,j,k,iblock))*WT1, &
                              ttt(i,j,k-1)*vit(i,j,k-1,iblock)+(1.0D0-vit(i,j,k-1,iblock))*WT1, &
                              ttt(i+1,j,k)*vit(i+1,j,k,iblock)+(1.0D0-vit(i+1,j,k,iblock))*WT1, &
                              ttt(i,j,k+1)*vit(i,j,k+1,iblock)+(1.0D0-vit(i,j,k+1,iblock))*WT1)
             atmin(i,j,k)=min(ttt(i,j,k)*vit(i,j,k,iblock)+(1.0D0-vit(i,j,k,iblock))*WT2, & 
                              ttt(i,j-1,k)*vit(i,j-1,k,iblock)+(1.0D0-vit(i,j-1,k,iblock))*WT2, & 
                              ttt(i,j+1,k)*vit(i,j+1,k,iblock)+(1.0D0-vit(i,j+1,k,iblock))*WT2, & 
                              ttt(i-1,j,k)*vit(i-1,j,k,iblock)+(1.0D0-vit(i-1,j,k,iblock))*WT2, &
                              ttt(i+1,j,k)*vit(i+1,j,k,iblock)+(1.0D0-vit(i+1,j,k,iblock))*WT2, &
                              ttt(i,j,k-1)*vit(i,j,k-1,iblock)+(1.0D0-vit(i,j,k-1,iblock))*WT2, &
                              ttt(i,j,k+1)*vit(i,j,k+1,iblock)+(1.0D0-vit(i,j,k+1,iblock))*WT2)
           end if
          ENDDO
        ENDDO
      ENDDO
   
      DO K = 1,km
        DO J = 2, jmt-1
          DO I = 2, imt-1
            if (at00(i,j,k)>atmax(i,j,k).or.at00(i,j,k)<atmin(i,j,k)) then
              adv_xy1(i,j,k)=-(ttt(i+1,j,k)-ttt(i,j,k))*        &
                                 abs(u_wface(i+1,j,k))*tarea_r(i,j,iblock)
              adv_xy2(i,j,k)= (ttt(i,j,k)-ttt(i-1,j,k))*        &
                                 abs(u_wface(i,j,k))*tarea_r(i,j,iblock)
              adv_xy3(i,j,k)=-(ttt(i,j+1,k)-ttt(i,j,k))*        &
                                 abs(v_sface(i,j,k))*tarea_r(i,j,iblock)
              adv_xy4(i,j,k)= (ttt(i,j,k)-ttt(i,j-1,k))*        &
                                 abs(v_sface(i,j-1,k))*tarea_r(i,j,iblock)
              adv_xy2(i+1,j,k)= (ttt(i+1,j,k)-ttt(i,j,k))*        &
                                 abs(u_wface(i+1,j,k))*tarea_r(i+1,j,iblock)
              adv_xy1(i-1,j,k)=-(ttt(i,j,k)-ttt(i-1,j,k))*        &
                                 abs(u_wface(i,j,k))*tarea_r(i-1,j,iblock)
              adv_xy4(i,j+1,k)= (ttt(i,j+1,k)-ttt(i,j,k))*        &
                                 abs(v_sface(i,j,k))*tarea_r(i,j+1,iblock)
              adv_xy3(i,j-1,k)=-(ttt(i,j,k)-ttt(i,j-1,k))*        &
                                 abs(v_sface(i,j-1,k))*tarea_r(i,j-1,iblock)
              if ( k == 1 ) then
                 adv_zb2(i,j,k)= 0.5*abs(WWW(I,J,k+1))*ODZP(k)*(ttt(I,J,k)-ttt(I,J,k+1))
                 adv_zb1(i,j,k+1)=-0.5*abs(WWW(I,J,k+1))*ODZP(k+1)*(ttt(I,J,k)-ttt(I,J,k+1))
              else if ( k==km) then
                 adv_zb2(i,j,k-1)= 0.5*abs(WWW(I,J,k))*ODZP(k-1)*(ttt(I,J,k-1)-ttt(I,J,k))
                 adv_zb1(i,j,k)=-0.5*abs(WWW(I,J,k  ))*ODZP(k)*(ttt(I,J,k-1)-ttt(I,J,k))
              else
                 adv_zb2(i,j,k-1)= 0.5*abs(WWW(I,J,k))*ODZP(k-1)*(ttt(I,J,k-1)-ttt(I,J,k))
                 adv_zb1(i,j,k)=-0.5*abs(WWW(I,J,k  ))*ODZP(k)*(ttt(I,J,k-1)-ttt(I,J,k))
                 adv_zb2(i,j,k)= 0.5*abs(WWW(I,J,k+1))*ODZP(k)*(ttt(I,J,k)-ttt(I,J,k+1))
                 adv_zb1(i,j,k+1)=-0.5*abs(WWW(I,J,k+1))*ODZP(k+1)*(ttt(I,J,k)-ttt(I,J,k+1))
              end if
            endif
          END DO
        END DO
      END DO
   
      adv_xx=-(adv_x0+adv_xy1+adv_xy2+adv_c1)
      adv_yy=-(adv_y0+adv_xy3+adv_xy4+adv_c2)
      adv_zz=-(adv_za+adv_zb1+adv_zb2+adv_zc)
      adv_tt= adv_xx+adv_yy+adv_zz
!LPF20160829
      ax(:,:,:,mtracer,iblock)= adv_xx(:,:,:) !/dble(NSS)
      ay(:,:,:,mtracer,iblock)= adv_yy(:,:,:) !/dble(NSS)
      az(:,:,:,mtracer,iblock)= adv_zz(:,:,:) !/dble(NSS)
      !ax(:,:,:,mtracer,iblock)= ax(:,:,:,mtracer,iblock) + adv_xx(:,:,:) !/dble(NSS)
      !ay(:,:,:,mtracer,iblock)= ay(:,:,:,mtracer,iblock) + adv_yy(:,:,:) !/dble(NSS)
      !az(:,:,:,mtracer,iblock)= az(:,:,:,mtracer,iblock) + adv_zz(:,:,:) !/dble(NSS)
      !ax(:,:,:,mtracer,iblock)= ax(:,:,:,mtracer,iblock) + adv_xx(:,:,:)/dble(NSS)
      !ay(:,:,:,mtracer,iblock)= ay(:,:,:,mtracer,iblock) + adv_yy(:,:,:)/dble(NSS)
      !az(:,:,:,mtracer,iblock)= az(:,:,:,mtracer,iblock) + adv_zz(:,:,:)/dble(NSS)
!LPF20160829
!
      deallocate ( adv_c1,adv_c2,adv_zc)
      deallocate ( adv_x0, adv_y0, adv_xx, adv_yy)
      deallocate ( adv_xy1, adv_xy2, adv_xy3, adv_xy4)
      deallocate ( at00, atmax, atmin )
      deallocate (adv_zz, adv_za, adv_zb1 )
      deallocate (adv_zb2, atmaxz, atminz, atz)
!
      else
        call exit_licom(sigAbort,'The false advection option for tracer')
      end if
!

      end subroutine advection_tracer

!***********************************************************************

 end module advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
