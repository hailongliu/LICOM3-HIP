!!     =================
!      SUBROUTINE JRA_DAILY(TNUM)
      SUBROUTINE JRA_DAILY
!     =================

#include <def-undef.h>
 use precision_mod
 use param_mod
 use pconst_mod
 use forc_mod
 use constant_mod
 use grid
 use domain
 use gather_scatter
 use POP_GridHorzMod
 use POP_HaloMod
! use dyn_mod, only: u,v,buffer,buffer_s
 use dyn_mod, only: u,v,buffer
 use tracer_mod, only: at
 use msg_mod
 use blocks

      IMPLICIT NONE
#include <netcdf.inc>

      integer :: mon_day,irec,tnum
      integer ::  ErrorCode
!      real(r8),dimension(s_imt,s_jmt,365) :: t10,u10,v10,slp,q10,swhf,lwhf
!      real(r8),dimension(s_imt,s_jmt,365) :: precr,precs,rf,si
! avoid > 2GB common for 5km version
!      real(r8),allocatable,dimension(:,:,:) :: t10,u10,v10,slp,q10,swhf,lwhf
!      real(r8),allocatable,dimension(:,:,:) :: precr,precs,rf,si

      real(r8),dimension(imt,jmt) :: model_sst,es,qs,zz,uu,vv,windx,windy,theta
      real(r8),dimension(imt,jmt) :: core_sensible,core_latent,core_tau
!
      real(r8), parameter :: tok=273.15
      real(r8), parameter :: epsln=1e-25
      real(r8),dimension(imt,jmt) :: tmp1,tmp2
      real(r8),dimension(imt,jmt,1) :: tmp1_su,tmp2_sv
      real(r8),dimension(imt_global,jmt_global) :: tmp3
      integer,save:: jra_init=0

!      if ( TNUM.gt.1 ) goto 111

! decide recode number
!      IYFM,MON0,IDAY
      mon_day=0
      do i=1,mon0-1
      mon_day=mon_day+nmonth(i)
      enddo
!
! start from a 1000-year spinup
!      irec=(iyfm-1)*365+mon_day+iday-365*(13-1)
!climatology forcing
      irec=mon_day+iday

      if (mytid.eq.0) then
!      write(*,*) "iyfm=",iyfm
!      write(*,*) "mon0=",mon0
!      write(*,*) "iday=",iday
!      write(*,*) "mon_day=",mon_day
      write(*,*) "irec=",irec
!      write(*,*) s_imt,s_jmt
      endif

! read in jra data
! note the dimensions are s_imt, s_jmt
!#ifdef SPMD
      if (jra_init.eq.0) then
      jra_init=1

      if (mytid.eq.master_task) then
      allocate(t10(s_imt,s_jmt,365))
      allocate(u10(s_imt,s_jmt,365))
      allocate(v10(s_imt,s_jmt,365))
      allocate(slp(s_imt,s_jmt,365))
      allocate(q10(s_imt,s_jmt,365))
      allocate(swhf(s_imt,s_jmt,365))
      allocate(lwhf(s_imt,s_jmt,365))
      allocate(precr(s_imt,s_jmt,365))
      allocate(precs(s_imt,s_jmt,365))
      allocate(rf(s_imt,s_jmt,365))
      allocate(si(s_imt,s_jmt,365))
      call read_jra(irec,"t_10.2016_daily.nc",t10)
      call read_jra(irec,"u_10.2016_daily.nc",u10)
      call read_jra(irec,"v_10.2016_daily.nc",v10)
      call read_jra(irec,"slp.2016_daily.nc",slp)
      call read_jra(irec,"q_10.2016_daily.nc",q10)
      call read_jra(irec,"rsds.2016_daily.nc",swhf)
      call read_jra(irec,"rlds.2016_daily.nc",lwhf)
      call read_jra(irec,"rain.2016_daily.nc",precr)
      call read_jra(irec,"snow.2016_daily.nc",precs)
      call read_jra(irec,"runoff_all.2016_daily.nc",rf)
      call read_jra1(irec,"ice.2016_daily.nc",si)

      allocate(buffer3d(imt,jmt,km,nblocks_tot))
!      allocate(w3d(imt_global,jmt_global,km))

!      print*,s_lon(:,1)
!      print*,s_lat(1,:)
!      stop 888888
!      print*,lon_o(:,1)
!      print*,lat_o(1,:)
!      stop 888888
!       print*,si
!       stop
      endif
      if(mytid==master_task+1) allocate(t10(s_imt,s_jmt,365))
      if(mytid==master_task+2) allocate(u10(s_imt,s_jmt,365))
      if(mytid==master_task+3) allocate(v10(s_imt,s_jmt,365))
      if(mytid==master_task+4) allocate(slp(s_imt,s_jmt,365))
      if(mytid==master_task+5) allocate(q10(s_imt,s_jmt,365))
      if(mytid==master_task+6) allocate(swhf(s_imt,s_jmt,365))
      if(mytid==master_task+7) allocate(lwhf(s_imt,s_jmt,365))
      if(mytid==master_task+8) allocate(precr(s_imt,s_jmt,365))
      if(mytid==master_task+9) allocate(precs(s_imt,s_jmt,365))
      if(mytid==master_task+10) allocate(rf(s_imt,s_jmt,365))
      if(mytid==master_task+11) allocate(si(s_imt,s_jmt,365))

      if(mytid==master_task+50) allocate(buffer3d(imt,jmt,km,nblocks_tot))
      if(mytid==master_task+20) allocate(buffer3d(imt,jmt,km,nblocks_tot))
      if(mytid==master_task+30) allocate(buffer3d(imt,jmt,km,nblocks_tot))
      if(mytid==master_task+40) allocate(buffer3d(imt,jmt,km,nblocks_tot))

      call MPI_Bcast(s_lon,s_imt*s_jmt,MPI_PR,master_task,mpi_comm_ocn,ierr)
      call MPI_Bcast(s_lat,s_imt*s_jmt,MPI_PR,master_task,mpi_comm_ocn,ierr)
      call MPI_Bcast(lon_o,imt_global*jmt_global,MPI_real4,master_task,mpi_comm_ocn,ierr)
      call MPI_Bcast(lat_o,imt_global*jmt_global,MPI_real4,master_task,mpi_comm_ocn,ierr)
      endif

!      call MPI_Bcast(t10(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,mpi_comm_ocn,ierr)
!      call MPI_Bcast(u10(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,mpi_comm_ocn,ierr)
!      call MPI_Bcast(v10(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,mpi_comm_ocn,ierr)
!      call MPI_Bcast(slp(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,mpi_comm_ocn,ierr)
!      call MPI_Bcast(q10(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,mpi_comm_ocn,ierr)
!      call MPI_Bcast(swhf(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,mpi_comm_ocn,ierr)
!      call MPI_Bcast(lwhf(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,mpi_comm_ocn,ierr)
!      call MPI_Bcast(precr(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,mpi_comm_ocn,ierr)
!      call MPI_Bcast(precs(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,mpi_comm_ocn,ierr)
!      call MPI_Bcast(rf(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,mpi_comm_ocn,ierr)
!      call MPI_Bcast(si(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,mpi_comm_ocn,ierr)

      if(mytid==master_task) then
      call MPI_send(t10(:,:,irec),s_imt*s_jmt,MPI_PR,master_task+1,123456,mpi_comm_ocn,ierr)
      call MPI_send(u10(:,:,irec),s_imt*s_jmt,MPI_PR,master_task+2,123456,mpi_comm_ocn,ierr)
      call MPI_send(v10(:,:,irec),s_imt*s_jmt,MPI_PR,master_task+3,123456,mpi_comm_ocn,ierr)
      call MPI_send(slp(:,:,irec),s_imt*s_jmt,MPI_PR,master_task+4,123456,mpi_comm_ocn,ierr)
      call MPI_send(q10(:,:,irec),s_imt*s_jmt,MPI_PR,master_task+5,123456,mpi_comm_ocn,ierr)
      call MPI_send(swhf(:,:,irec),s_imt*s_jmt,MPI_PR,master_task+6,123456,mpi_comm_ocn,ierr)
      call MPI_send(lwhf(:,:,irec),s_imt*s_jmt,MPI_PR,master_task+7,123456,mpi_comm_ocn,ierr)
      call MPI_send(precr(:,:,irec),s_imt*s_jmt,MPI_PR,master_task+8,123456,mpi_comm_ocn,ierr)
      call MPI_send(precs(:,:,irec),s_imt*s_jmt,MPI_PR,master_task+9,123456,mpi_comm_ocn,ierr)
      call MPI_send(rf(:,:,irec),s_imt*s_jmt,MPI_PR,master_task+10,123456,mpi_comm_ocn,ierr)
      call MPI_send(si(:,:,irec),s_imt*s_jmt,MPI_PR,master_task+11,123456,mpi_comm_ocn,ierr)
      end if

       
      allocate (buffer(imt_global,jmt_global))
!      allocate (buffer_s(s_imt,s_jmt))

! interplate to T grid
!first parallel interp
      if (mytid.eq.master_task+1) then
      call MPI_recv(t10(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,123456,mpi_comm_ocn,status,ierr)
!      if (mytid.eq.master_task) then
!      call interplation_nearest(t10(:,:,irec),tmp3)
      call interplation_nearest(t10(:,:,irec),tmp3)
!      call interplation_nearest(buffer_s,tmp3)
       buffer=tmp3 !LPF20200321
      endif

      if (mytid.eq.master_task+2) then
      call MPI_recv(u10(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,123456,mpi_comm_ocn,status,ierr)
!      if (mytid.eq.master_task) then
      call interplation_nearest(u10(:,:,irec),tmp3)
!      call interplation_nearest(buffer_s,tmp3)
       buffer=tmp3 !LPF20200321
      endif

      if (mytid.eq.master_task+3) then
      call MPI_recv(v10(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,123456,mpi_comm_ocn,status,ierr)
!      if (mytid.eq.master_task) then
      call interplation_nearest(v10(:,:,irec),tmp3)
!      call interplation_nearest(buffer_s,tmp3)
       buffer=tmp3 !LPF20200321
      endif

      if (mytid.eq.master_task+4) then
      call MPI_recv(slp(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,123456,mpi_comm_ocn,status,ierr)
!      if (mytid.eq.master_task) then
      call interplation_nearest(slp(:,:,irec),tmp3)
!      call interplation_nearest(buffer_s,tmp3)
       buffer=tmp3 !LPF20200321
      endif

      if (mytid.eq.master_task+5) then
      call MPI_recv(q10(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,123456,mpi_comm_ocn,status,ierr)
!      if (mytid.eq.master_task) then
      call interplation_nearest(q10(:,:,irec),tmp3)
!      call interplation_nearest(buffer_s,tmp3)
       buffer=tmp3 !LPF20200321
      endif

      if (mytid.eq.master_task+6) then
      call MPI_recv(swhf(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,123456,mpi_comm_ocn,status,ierr)
!      if (mytid.eq.master_task) then
      call interplation_nearest(swhf(:,:,irec),tmp3)
!      call interplation_nearest(buffer_s,tmp3)
       buffer=tmp3 !LPF20200321
      endif

      if (mytid.eq.master_task+7) then
      call MPI_recv(lwhf(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,123456,mpi_comm_ocn,status,ierr)
!      if (mytid.eq.master_task) then
      call interplation_nearest(lwhf(:,:,irec),tmp3)
!      call interplation_nearest(buffer_s,tmp3)
       buffer=tmp3 !LPF20200321
      endif

      if (mytid.eq.master_task+8) then
      call MPI_recv(precr(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,123456,mpi_comm_ocn,status,ierr)
!      if (mytid.eq.master_task) then
      call interplation_nearest(precr(:,:,irec),tmp3)
!      call interplation_nearest(buffer_s,tmp3)
       buffer=tmp3 !LPF20200321
      endif

      if (mytid.eq.master_task+9) then
      call MPI_recv(precs(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,123456,mpi_comm_ocn,status,ierr)
!      if (mytid.eq.master_task) then
      call interplation_nearest(precs(:,:,irec),tmp3)
!      call interplation_nearest(buffer_s,tmp3)
       buffer=tmp3 !LPF20200321
      endif

      if (mytid.eq.master_task+10) then
      call MPI_recv(rf(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,123456,mpi_comm_ocn,status,ierr)
!      if (mytid.eq.master_task) then
      call interplation_nearest(rf(:,:,irec),tmp3)
!      call interplation_nearest(buffer_s,tmp3)
       buffer=tmp3 !LPF20200321
      endif

      if (mytid.eq.master_task+11) then
      call MPI_recv(si(:,:,irec),s_imt*s_jmt,MPI_PR,master_task,123456,mpi_comm_ocn,status,ierr)
!      if (mytid.eq.master_task) then
      call interplation_nearest(si(:,:,irec),tmp3)
!      call interplation_nearest(buffer_s,tmp3)
       buffer=tmp3 !LPF20200321
      endif

!second scatter
      call scatter_global(tsa3(:,:,1,:),buffer, master_task+1, distrb_clinic, &
!      call scatter_global(tsa3(:,:,1,:),buffer, master_task, distrb_clinic, &
                       field_loc_center, field_type_scalar)
      call POP_HaloUpdate(tsa3 , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

      call scatter_global(wspdu3(:,:,1,:),buffer, master_task+2, distrb_clinic, &
!      call scatter_global(wspdu3(:,:,1,:),buffer, master_task, distrb_clinic, &
                       field_loc_center, field_type_scalar)
      call POP_HaloUpdate(wspdu3 , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

      call scatter_global(wspdv3(:,:,1,:),buffer, master_task+3, distrb_clinic, &
!      call scatter_global(wspdv3(:,:,1,:),buffer, master_task, distrb_clinic, &
                       field_loc_center, field_type_scalar)
      call POP_HaloUpdate(wspdv3 , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

      call scatter_global(psa3(:,:,1,:),buffer, master_task+4, distrb_clinic, &
!      call scatter_global(psa3(:,:,1,:),buffer, master_task, distrb_clinic, &
                       field_loc_center, field_type_scalar)
      call POP_HaloUpdate(psa3 , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

      call scatter_global(qar3(:,:,1,:),buffer, master_task+5, distrb_clinic, &
!      call scatter_global(qar3(:,:,1,:),buffer, master_task, distrb_clinic, &
                       field_loc_center, field_type_scalar)
      call POP_HaloUpdate(qar3 , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

      call scatter_global(swv3(:,:,1,:),buffer, master_task+6, distrb_clinic, &
!      call scatter_global(swv3(:,:,1,:),buffer, master_task, distrb_clinic, &
                      field_loc_center, field_type_scalar)
      call POP_HaloUpdate(swv3 , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

      call scatter_global(lwv3(:,:,1,:),buffer, master_task+7, distrb_clinic, &
!      call scatter_global(lwv3(:,:,1,:),buffer, master_task, distrb_clinic, &
                       field_loc_center, field_type_scalar)
      call POP_HaloUpdate(lwv3 , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

      call scatter_global(rain3(:,:,1,:),buffer, master_task+8, distrb_clinic, &
!      call scatter_global(rain3(:,:,1,:),buffer, master_task, distrb_clinic, &
                       field_loc_center, field_type_scalar)
      call POP_HaloUpdate(rain3 , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

      call scatter_global(snow3(:,:,1,:),buffer, master_task+9, distrb_clinic, &
!      call scatter_global(snow3(:,:,1,:),buffer, master_task, distrb_clinic, &
                       field_loc_center, field_type_scalar)
      call POP_HaloUpdate(snow3 , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

      call scatter_global(runoff3(:,:,1,:),buffer, master_task+10, distrb_clinic, &
!      call scatter_global(runoff3(:,:,1,:),buffer, master_task, distrb_clinic, &
                       field_loc_center, field_type_scalar)
      call POP_HaloUpdate(runoff3 , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

      call scatter_global(seaice3(:,:,1,:),buffer, master_task+11, distrb_clinic, &
!      call scatter_global(seaice3(:,:,1,:),buffer, master_task, distrb_clinic, &
                       field_loc_center, field_type_scalar)
      call POP_HaloUpdate(seaice3 , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

!111    continue

         do j = jsm,jem
            do i = 2,imm
              uu(i,j)=vit(i,j,1,1)*(u(i,j,1,1)+u(i+1,j,1,1)+u(i,j-1,1,1)+u(i+1,j-1,1,1)) &
              /(viv(i,j,1,1)+viv(i+1,j,1,1)+viv(i,j-1,1,1)+viv(i+1,j-1,1,1)+epsln)
              vv(i,j)=vit(i,j,1,1)*(v(i,j,1,1)+v(i+1,j,1,1)+v(i,j-1,1,1)+v(i+1,j-1,1,1)) &
              /(viv(i,j,1,1)+viv(i+1,j,1,1)+viv(i,j-1,1,1)+viv(i+1,j-1,1,1)+epsln)
            end do
         end do

! transfer core data to what the subroutine need
      do j=1,jmt
         do i=1,imt
! relative speed to surface currents
         windx(i,j)=(wspdu3(i,j,1,1)-uu(i,j))*vit(i,j,1,1)
         windy(i,j)=(wspdv3(i,j,1,1)+vv(i,j))*vit(i,j,1,1)
!         windy(i,j)=(wspdv3(i,j,1,1)-vv(i,j))*vit(i,j,1,1)
!  1.0 is from mom4
         wspd3(i,j,1,1)=sqrt(windx(i,j)**2+windy(i,j)**2+1.0)*vit(i,j,1,1)
! using a transient temperature, not daily mean
         model_sst(i,j)=(at(i,j,1,1,1)+tok)*vit(i,j,1,1)
         zz(i,j)=10.
         qs(i,j)=0.98*640380*exp(-5107.4/model_sst(i,j))/1.22*vit(i,j,1,1)
! temperature to potential temperature
         theta(i,j)=tsa3(i,j,1,1)*(100000.0/psa3(i,j,1,1))**0.286*vit(i,j,1,1)
         runoff(i,j,1)=runoff3(i,j,1,1)*vit(i,j,1,1)
         seaice(i,j,1)=seaice3(i,j,1,1)*vit(i,j,1,1)
         end do
      end do


! compute heat flux
       call ncar_ocean_fluxes_jra(wspd3(1,1,1,1),theta(1,1),model_sst(1,1),qar3(1,1,1,1),qs(1,1),zz(1,1),vit(1,1,1,1),&
           core_sensible(1,1),core_latent(1,1),core_tau,ustar(1,1,1))

       do j=1,jmt
          do i=1,imt
            sshf(i,j,1)=core_sensible(i,j)*vit(i,j,1,1)*(1.0d0-seaice(i,j,1))
            lthf(i,j,1)=(core_latent(i,j)-snow3(i,j,1,1)*3.335e+5)*vit(i,j,1,1)*(1.0d0-seaice(i,j,1))
            lwv(i,j,1)=(0.95*lwv3(i,j,1,1)-0.95*5.67E-8*model_sst(i,j)**4)*vit(i,j,1,1)*(1.0d0-seaice(i,j,1))
            tmp1(i,j)=core_tau(i,j)*windx(i,j)*vit(i,j,1,1)
            tmp2(i,j)=-core_tau(i,j)*windy(i,j)*vit(i,j,1,1)
            nswv(i,j,1)=(lwv(i,j,1)+sshf(i,j,1)+lthf(i,j,1))
            swv(i,j,1)=(1-0.066)*swv3(i,j,1,1)*vit(i,j,1,1)*(1.0d0-seaice(i,j,1))
            fresh(i,j,1)=-(core_latent(i,j)/(2.5e+6)+rain3(i,j,1,1)+snow3(i,j,1,1)+runoff(i,j,1))*vit(i,j,1,1)*(1.0d0-seaice(i,j,1))
            ustar(i,j,1)=ustar(i,j,1)*vit(i,j,1,1)
         end do
      end do
!multiply the anglet !updated 20200520
        do j = jst,jem
           do i = 2,imm
           tmp1_su(i,j,1)= tmp1(i,j)*cos(anglet(i,j,1))-tmp2(i,j)*sin(anglet(i,j,1))
           !tmp1_su(i,j,1)= tmp1(i,j)*cos(anglet(i,j,1))+tmp2(i,j)*sin(anglet(i,j,1))
           !tmp2_sv(i,j,1)=-tmp2(i,j)*cos(anglet(i,j,1))+tmp1(i,j)*sin(anglet(i,j,1))
           tmp2_sv(i,j,1)=tmp2(i,j)*cos(anglet(i,j,1))+tmp1(i,j)*sin(anglet(i,j,1))
           end do
        end do
!update the halo ??
 
! tau to U/V grid
        do j = jst,jem
           do i = 2,imm
             su(i,j,1)= 0.25*(tmp1_su(i,j,1)+tmp1_su(i-1,j,1)+tmp1_su(i,j+1,1)+tmp1_su(i-1,j+1,1))*viv(i,j,1,1)
             sv(i,j,1)= 0.25*(tmp2_sv(i,j,1)+tmp2_sv(i-1,j,1)+tmp2_sv(i,j+1,1)+tmp2_sv(i-1,j+1,1))*viv(i,j,1,1)
             !su(i,j,1)= 0.25*(tmp1(i,j)+tmp1(i-1,j)+tmp1(i,j+1)+tmp1(i-1,j+1))*viv(i,j,1,1)
             !sv(i,j,1)= 0.25*(tmp2(i,j)+tmp2(i-1,j)+tmp2(i,j+1)+tmp2(i-1,j+1))*viv(i,j,1,1)
           !tmp_su(i,j,iblock)= taux(i,j,iblock)*cos(anglet(i,j,iblock))+tauy(i,j,iblock)*sin(anglet(i,j,iblock))
           !tmp_sv(i,j,iblock)=-tauy(i,j,iblock)*cos(anglet(i,j,iblock))+taux(i,j,iblock)*sin(anglet(i,j,iblock))
           end do
        end do

      deallocate(buffer)
!      deallocate(buffer_s)
      return
      end subroutine JRA_DAILY

!---------------------------------------------
      subroutine read_jra(nnn,fname,var)
!---------------------------------------------
use precision_mod
      use param_mod, only: s_imt,s_jmt
      use pconst_mod, only: s_lon,s_lat
      implicit none
#include <netcdf.inc>

      integer :: start(3),count(3)
      integer :: ncid,iret,nnn,i,j
      real(r8) :: var(s_imt,s_jmt,365)
      real(r4) :: ttmp(s_imt,s_jmt,365)
      character (len=180) :: fname

      start(1)=1;count(1)=s_imt
      start(2)=1;count(2)=s_jmt
!      start(3)=nnn;count(3)=1
      start(3)=1;count(3)=365

      iret=nf_open(fname,nf_nowrite,ncid)
      call check_err (iret)
      iret=nf_get_var_double(ncid,3,s_lon(:,1))
      call check_err (iret)
      iret=nf_get_var_double(ncid,4,s_lat(1,:))
      call check_err (iret)
      iret=nf_get_vara_real(ncid,5,start,count,ttmp)
      call check_err (iret)

      iret=nf_close(ncid)
      call check_err (iret)

      do j=2,s_jmt
      do i=1,s_imt
      s_lon(i,j)=s_lon(i,1)
      enddo
      enddo
      do j=1,s_jmt
      do i=2,s_imt
      s_lat(i,j)=s_lat(1,j)
      enddo
      enddo
!      print*,s_lon(:,1)
!      print*,s_lat(1,:)
!      stop

      var=ttmp

      return
      end subroutine read_jra

!---------------------------------------------
      subroutine read_jra1(nnn,fname,var)
!---------------------------------------------
use precision_mod
      use param_mod, only: s_imt,s_jmt
      use pconst_mod, only: s_lon,s_lat
      implicit none
#include <netcdf.inc>

      integer :: start(3),count(3)
      integer :: ncid,iret,nnn,i,j
      real(r8) :: var(s_imt,s_jmt,365)
      real(r4) :: ttmp(s_imt,s_jmt,365),xt(s_imt),yt(s_jmt)
      character (len=180) :: fname

      start(1)=1;count(1)=s_imt
      start(2)=1;count(2)=s_jmt
!      start(3)=nnn;count(3)=1
      start(3)=1;count(3)=365

      iret=nf_open(fname,nf_nowrite,ncid)
      call check_err (iret)
      iret=nf_get_var_real(ncid,3,xt)
      call check_err (iret)
      iret=nf_get_var_real(ncid,4,yt)
      call check_err (iret)
      iret=nf_get_vara_real(ncid,5,start,count,ttmp)
      call check_err (iret)

      iret=nf_close(ncid)
      call check_err (iret)

      do j=2,s_jmt
      do i=1,s_imt
      s_lon(i,j)=xt(i)
      enddo
      enddo
      do j=1,s_jmt
      do i=2,s_imt
      s_lat(i,j)=yt(j)
      enddo
      enddo
!      print*,s_lon(:,1)
!      print*,s_lat(1,:)
!      stop

      var=ttmp

      return
      end subroutine read_jra1

!---------------------------------------------
      subroutine interplation_nearest(source,object)
!---------------------------------------------
use precision_mod
      use param_mod, only: imt,jmt,imt_global,jmt_global,s_imt,s_jmt,mytid
      use pconst_mod, only: s_lon,s_lat
      use output_mod, only: spval
      implicit none

      integer, parameter :: iwk=s_imt+2,jwk=s_jmt+2
      real(r8), parameter :: dx=0.5625,dy=0.5625
      real(r8) :: source(s_imt,s_jmt)
      real(r8) :: object(imt_global,jmt_global)
      real(r8) :: s_work(iwk,jwk),s_wx(iwk,jwk),s_wy(iwk,jwk)
      integer :: i,j

      do j=1,s_jmt
      do i=1,s_imt
      s_wx(i+1,j+1)=s_lon(i,j)
      enddo
      enddo

      do j=1,s_jmt
      s_wx(  1,j+1)=s_wx(    2,j+1)-dx
      s_wx(iwk,j+1)=s_wx(iwk-1,j+1)+dx
      enddo
!
      do i=1,s_imt
      s_wx(i+1,  1)=s_wx(i+1,2)
      s_wx(i+1,jwk)=s_wx(i+1,1)
      enddo

      s_wx(  1,  1)=s_wx(    1,    2)
      s_wx(iwk,  1)=s_wx(  iwk,    2)
      s_wx(  1,jwk)=s_wx(    1,jwk-1)
      s_wx(iwk,jwk)=s_wx(  iwk,jwk-1)

      do j=1,s_jmt
      do i=1,s_imt
      s_wy(i+1,j+1)=s_lat(i,j)
      enddo
      enddo

      do j=1,s_jmt
      s_wy(  1,j+1)=s_wy(    2,j+1)
      s_wy(iwk,j+1)=s_wy(iwk-1,j+1)
      enddo
      s_wy(  1,  1)=s_wy(  1,  2)+dy
      s_wy(iwk,  1)=s_wy(iwk,  2)+dy
!
      do i=1,s_imt
      s_wy(i+1,  1)=s_wy(i+1,    2)+dy
      s_wy(i+1,jwk)=s_wy(i+1,jwk-1)-dy
      enddo
      s_wy(  1,jwk)=s_wy(    1,jwk-1)-dy
      s_wy(iwk,jwk)=s_wy(  iwk,jwk-1)-dy
!

      do j=1,s_jmt
      do i=1,s_imt
      s_work(i+1,j+1)=source(i,j)
      enddo
      enddo

      do j=1,s_jmt
      s_work(  1,j+1)=s_work(    2,j+1)
      s_work(iwk,j+1)=s_work(iwk-1,j+1)
      enddo
!
      do i=1,iwk
      s_work(i,  1)=s_work(i,    2)
      s_work(i,jwk)=s_work(i,jwk-1)
      enddo
!
      call near(s_work,iwk,jwk,s_wx,s_wy,object(1,1))

      return
      end subroutine interplation_nearest

!---------------------------------------------
      subroutine near(a,mx,my,alon,alat,b)
!---------------------------------------------
!     input : a
!     output: b
!     A bi-linear interpolation will be used for the initial guess, then
!     a refilling procedure be called to redefine the missing data.
!     mx,my            x and y grids number of a
!     alon,alat        lontitude and latitude of a
!     spval            missing flag
!     nf               1 for T grid; 0 for U grid
!
use precision_mod
      use param_mod, only: imt_global,jmt_global,s_imt,s_jmt
      use pconst_mod, only: lon_o,lat_o!vit !,vit_global_surface
      use output_mod, only: spval
      implicit none
!
      integer :: mx,my,ic,jc,ip,jp,i,j
      real(r8):: b(imt_global,jmt_global)
      real(r8):: a(mx,my),alat(mx,my),alon(mx,my)
      real(r8), parameter :: isp=99999.0d0
      real(r8):: lon_tmp
!
!   Initiallize
!
!  ---nf = 1 on T girds

!  ---b results
      do j=1,jmt_global
      do i=1,imt_global
      b(i,j)=spval
      enddo
      enddo
!
      do 100 j=1,jmt_global
      do 100 i=1,imt_global
!
!lhl      if(vit(i,j,1,1).lt.0.01) goto 100
!
!  ---find adjacent two grids on x direction (ic)
      ic=isp
      jc=isp

      if (lon_o(i,j).lt.0.0) then
      lon_tmp=lon_o(i,j)+360.
      else
      lon_tmp=lon_o(i,j)
      endif

      do 45 ip=2,mx
      if(alon(ip-1,1).le.lon_tmp.and.alon(ip,1).ge.lon_tmp)then
      ic=ip
      goto 33
      endif
   45 continue
   33 continue

!
!  ---find adjacent two grids on y direction (jc)
      do 50 jp=2,my
      if(alat(1,jp).le.lat_o(i,j).and.alat(1,jp-1).ge.lat_o(i,j))then
      jc=jp
      goto 44
      endif
   50 continue
   44 continue

!
!  ---break if adjacent grids has no found
      if(ic.eq.isp.or.jc.eq.isp) then
      write(*,*)ic,jc
      write(*,*)i,j,lon_o(i,j),lat_o(i,j)
      write(*,*) alon(:,1)
      write(*,*) alat(1,:)
      stop 999
      endif
!
!  Bilinear interpolater
!
      b(i,j)=a(ic,jc)

!lhl      endif
!
  100 continue
!
!    set cyclic b. c.
!
      do j=1,jmt_global
      b(imt_global-1,j)     = b(1,j)
      b(imt_global,j)       = b(2,j)
      enddo
!
      return
      end subroutine near

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Over-ocean fluxes following Large and Yeager (used in NCAR models)           !
! Coded by Mike Winton (Michael.Winton@noaa.gov) in 2004
!
! A bug was found by Laurent Brodeau (brodeau@gmail.com) in 2007.
! Stephen.Griffies@noaa.gov updated the code with the bug fix. 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!
subroutine ncar_ocean_fluxes_jra (u_del, t, ts, q, qs, z, avail, &
                              sh,lh,tau,ustar)
!                              cd, ch, ce, ustar, bstar       )
use precision_mod
      use param_mod, only: imt,jmt,imt_global,jmt_global,s_imt,s_jmt,mytid
      implicit none
    real(r8)   , intent(in)   , dimension(imt,jmt) :: u_del, t, ts, q, qs, z
    real(r8)   , intent(in)   , dimension(imt,jmt) :: avail
    real(r8)   , intent(inout), dimension(imt,jmt) :: lh,sh,tau
    real(r8)   , dimension(imt,jmt) :: cd, ch, ce, ustar, bstar
!    real   , intent(inout), dimension(imt,jmt) :: cd, ch, ce, ustar, bstar

  real(r8) :: cd_n10, ce_n10, ch_n10, cd_n10_rt    ! neutral 10m drag coefficients
  real(r8) :: cd_rt                                ! full drag coefficients @ z
  real(r8) :: zeta, x2, x, psi_m, psi_h            ! stability parameters
  real(r8) :: u, u10, tv, tstar, qstar, z0, xx, stab

  integer, parameter :: n_itts = 2
  real(r8), parameter :: grav = 9.80, vonkarm = 0.40,  L=2.5e6, cp=1000.5, r0=1.22
  integer               i, j, jj


  do j=1,jmt
  do i=1,imt
!  do i=1,size(u_del)
    if (avail(i,j) > 0.5 ) then
      tv = t(i,j)*(1+0.608*q(i,j));
      u = max(u_del(i,j), 0.5);                                 ! 0.5 m/s floor on wind (undocumented NCAR)
      u10 = u;                                                ! first guess 10m wind
    
      cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                ! L-Y eqn. 6a
      cd_n10_rt = sqrt(cd_n10);
      ce_n10 =                     34.6 *cd_n10_rt/1e3;       ! L-Y eqn. 6b
      stab = 0.5 + sign(0.5,t(i,j)-ts(i,j))
      ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;       ! L-Y eqn. 6c
  
      cd(i,j) = cd_n10;                                         ! first guess for exchange coeff's at z
      ch(i,j) = ch_n10;
      ce(i,j) = ce_n10;
      do jj=1,n_itts                                           ! Monin-Obukhov iteration
        cd_rt = sqrt(cd(i,j));
        ustar(i,j) = cd_rt*u;                                   ! L-Y eqn. 7a
        tstar    = (ch(i,j)/cd_rt)*(t(i,j)-ts(i,j));                ! L-Y eqn. 7b
        qstar    = (ce(i,j)/cd_rt)*(q(i,j)-qs(i,j));                ! L-Y eqn. 7c
        bstar(i,j) = grav*(tstar/tv+qstar/(q(i,j)+1/0.608));

        zeta     = vonkarm*bstar(i,j)*z(i,j)/(ustar(i,j)*ustar(i,j)); ! L-Y eqn. 8a
        zeta     = sign( min(abs(zeta),10.0), zeta );         ! undocumented NCAR
        x2 = sqrt(abs(1-16*zeta));                            ! L-Y eqn. 8b
        x2 = max(x2, 1.0);                                    ! undocumented NCAR
        x = sqrt(x2);
    
        if (zeta > 0) then
          psi_m = -5*zeta;                                    ! L-Y eqn. 8c
          psi_h = -5*zeta;                                    ! L-Y eqn. 8c
        else
          psi_m = log((1+2*x+x2)*(1+x2)/8)-2*(atan(x)-atan(1.0)); ! L-Y eqn. 8d
          psi_h = 2*log((1+x2)/2);                                ! L-Y eqn. 8e
        end if
    
        u10 = u/(1+cd_n10_rt*(log(z(i,j)/10)-psi_m)/vonkarm);       ! L-Y eqn. 9


        cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                  ! L-Y eqn. 6a again
        cd_n10_rt = sqrt(cd_n10);
        ce_n10 = 34.6*cd_n10_rt/1e3;                              ! L-Y eqn. 6b again
        stab = 0.5 + sign(0.5,zeta)
        ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;         ! L-Y eqn. 6c again
        z0 = 10*exp(-vonkarm/cd_n10_rt);                          ! diagnostic
    
        xx = (log(z(i,j)/10)-psi_m)/vonkarm;
        cd(i,j) = cd_n10/(1+cd_n10_rt*xx)**2;                       ! L-Y 10a
        xx = (log(z(i,j)/10)-psi_h)/vonkarm;
!!$        ch(i,j) = ch_n10/(1+ch_n10*xx/cd_n10_rt)**2;                !       b (bug)
!!$        ce(i,j) = ce_n10/(1+ce_n10*xx/cd_n10_rt)**2;                !       c (bug)
        ch(i,j) = ch_n10/(1+ch_n10*xx/cd_n10_rt)*sqrt(cd(i,j)/cd_n10) ! 10b (corrected code aug2007)
        ce(i,j) = ce_n10/(1+ce_n10*xx/cd_n10_rt)*sqrt(cd(i,j)/cd_n10) ! 10c (corrected code aug2007)
      end do
    end if
  end do
  end do

    sh=r0*cp*ch*(t-ts)*u_del
    lh=r0*ce*l*(q-qs)*u_del
    tau=r0*cd*u_del

    return
end subroutine ncar_ocean_fluxes_jra
