#define LOGMSG()
!write(mytid+600,'(a,i4)')"LICOM",__LINE__
#define LOGMSGCarbon()
!write(600+mytid,'(a,i4)')"Licom CARBON",__LINE__

#ifdef COUP
   use mct_mod
   use shr_dmodel_mod
   use esmf
   use seq_flds_mod
   use seq_cdata_mod
   use seq_infodata_mod
   use seq_timemgr_mod
   use shr_file_mod
   use shr_cal_mod, only : shr_cal_date2ymd
   use perf_mod
   use shr_dmodel_mod
   use fluxcpl
   use POP_CplIndices
   use shr_msg_mod
   use shr_cal_mod,       only: shr_cal_date2ymd
#endif
   use shr_sys_mod
   use shr_const_mod,only:SHR_CONST_SPVAL !linpf 2012Jul26
   use POP_CommMod
   use domain
   use grid
   use blocks
   use gather_scatter
   use distribution

#include <def-undef.h>
use param_mod
use pconst_mod
use control_mod
!use control_licom_mod
use constant_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
use tracer_mod
use pmix_mod
use forc_mod
!
#ifdef USE_OCN_CARBON
use carbon_mod
use cforce_mod
#endif


  implicit none
#include <netcdf.inc>

  integer :: mm
  REAL(r8):: t3a,t3b,t1,t2,clock0f
external allocate_fortran
external allocate_common

  mpi_comm_ocn=0

  cdate = 0
  sec = 0
  ierr = 0
  info_time = 0
      call mpi_init (ierr)
      call mpi_comm_dup(mpi_comm_world, mpi_comm_ocn, ierr)
      call mpi_comm_rank (mpi_comm_ocn, mytid, ierr)
      if (mytid == 0 ) then
!        open (15,file="ocn_modelio.nml", form="formatted")
!        read (15,modelio)
!        close(15)
         open (16,file="ocn.log", form="formatted")
      end if
!
      if (mytid==0)  write(16,*)"End mpi_comm_rank"
      call mpi_comm_size (mpi_comm_ocn, nproc, ierr)
      if (mytid ==0) write(16,*) "MYTID=",mytid,"Number of Processors is",nproc
!YU

      my_task = mytid
      master_task = 0
      POP_mytask = mytid
      POP_mastertask = 0


!---------------------------------------------------------------------
!     SET THE CONSTANTS USED IN THE MODEL
!---------------------------------------------------------------------
#ifdef COUP
      call shr_sys_flush(6)
#endif
    LOGMSG()
      CALL CONST
      call allocate_common
    LOGMSG()
      if (mytid == 0) then
      write(111,*)"OK------3"
      close(111)
      end if
#ifdef COUP
      call shr_sys_flush(6)
#endif
#ifdef SHOW_TIME
      call run_time('CONST')
#endif

!***********************************************************************
!          SET SOME CONSTANTS FOR THE BOGCM
!***********************************************************************
#ifdef USE_OCN_CARBON
      LOGMSGCarbon()
      CALL CTRLC
      LOGMSGCarbon()
#ifdef SHOW_TIME
      call run_time('CTRLC')
#endif
#endif

!---------------------------------------------------------------------
!     SET MODEL'S RESOLUTION,TOPOGRAPHY AND THE CONSTANT
!     PARAMETERS RELATED TO LATITUDES (J)
!---------------------------------------------------------------------
    LOGMSG()
   call init_domain_blocks
      if (mytid == 0) then
      write(111,*)"OK------3"
      close(111)
      end if
   call init_grid1
      if (mytid == 0) then
      write(111,*)"OK------3.1"
      close(111)
      end if
   call init_domain_distribution(KMT_G)
      if (mytid == 0) then
      write(111,*)"OK------3.2"
      close(111)
      end if
   call init_grid2
      if (mytid == 0) then
      write(111,*)"OK------3.3"
      close(111)
      end if
!      if(1) call SSAVELATLON !LPF20200313

   CALL GRIDS
      if (mytid == 0) then
      write(111,*)"OK------3.4"
      close(111)
      end if
   call calc_coeff
      if (mytid == 0) then
      write(111,*)"OK------4"
      close(111)
      end if
#ifdef SHOW_TIME
      call run_time('GRIDS')
#endif


!---------------------------------------------------------------------
!     SET SURFACE FORCING FIELDS (1: Annual mean; 0: Seasonal cycle)
!---------------------------------------------------------------------
    LOGMSG()
     CALL RDRIVER !LPF20151006
#ifdef SHOW_TIME
      call run_time('RDRIVER')
#endif
      if (mytid == 0) then
      write(111,*)"OK------5"
      close(111)
      end if

!***********************************************************************
!      SET FORCING DATA USED IN CARBON CYCLYE
!***********************************************************************
#ifdef USE_OCN_CARBON
      LOGMSGCarbon()
       CALL CFORCE
      LOGMSGCarbon()
#ifdef SHOW_TIME
       call run_time('CFORCE')
#endif
#endif

!---------------------------------------------------------------------
!     INITIALIZATION
!---------------------------------------------------------------------
    LOGMSG()
      CALL INIRUN
      if (mytid == 0) then
      write(111,*)"OK------5.0"
      close(111)
      end if
#ifdef SHOW_TIME
      call run_time('INRUN')
#endif
    if (NSTART==1)then !initial run from observation T/S
      number_month=1
      number_day=1
      yy_licom=int(number_month/12)+1
      mm_licom=number_month-(yy_licom-1)*12
      dd_licom=number_day
      curr_ymd_licom=10000*yy_licom+100*mm_licom+dd_licom
      iyfm = 1
      month = number_month

!   IF (mytid.eq.0) THEN
!     if (curr_ymd_licom.ne.curr_ymd_cpl) then
!      if (mytid==0) write (*,*) "initial state mismacth! STOP Now"
!         stop 
!     endif
!   ENDIF
   else
!LPF20131027 adjusted for restart 12
      if (mod(number_month,12)==0) then
         yy_licom=int(number_month/12)
         mm_licom=number_month-(yy_licom-1)*12
         dd_licom=number_day
!LPF20131027 adjusted for restart 12
      else
         yy_licom=int(number_month/12)+1
         mm_licom=number_month-(yy_licom-1)*12
         dd_licom=number_day
     endif
   endif
      month = mm_licom  !number_month
      iyfm=yy_licom
      if (mytid==0) write (*,*) "month,number_month,mm_licom",month,number_month,mm_licom
!lhl20130419


#ifdef USE_OCN_CARBON
      LOGMSGCarbon()
      CALL INIRUN_PT
      LOGMSGCarbon()
#ifdef SHOW_TIME
      call run_time('INIRUN_PT')
#endif
#endif


#ifdef CANUTO
      call turb_ini
#endif
!---------------------------------------------------------------------
!     INITIALIZATION OF ISOPYCNAL MIXING
!---------------------------------------------------------------------
    LOGMSG()
#ifdef ISO
      CALL ISOPYI
#ifdef SHOW_TIME
      call run_time('ISOPYI')
#endif
#endif

#ifdef COUP
         call shr_sys_flush(6)
#endif

!linpf 2012Jul27
!----------------------------------------------------------------------
!   set time 
!----------------------------------------------------------------------
#if (defined SOLARCHLORO) 
!!$OMP PARALLEL DO PRIVATE (J,I)
        DO J = JST,JET
         DO I = 1,IMT
           chloro(I,J)= ( (chloro3 (I,J,IPT2) - chloro3 (I,J,IPT1))     &
                    * FACTOR + chloro3 (I,J,IPT1))
           END DO
        END DO

          CALL  SW_ABSOR
#endif
!LPF20131024

!       if (mytid.eq.0) print*,sss(1,:)
!       if (mytid.eq.0) print*,sss3(1,:,IPT1)
!       if (mytid.eq.0) print*,sss3(1,:,IPT2)
           
!---------------------------------------------------------------------
!     THERMAL CYCLE
!---------------------------------------------------------------------
      CALL allocate_fortran
       if(mytid==0) write(*,*) 'stepon begin, version: 05km (GPU/HIP)'
#define GPU
#ifdef COUP
       DO II = 1,NSS/num_cpl !thermal cycle in a couple interval
#else
       DO MM = 1, number
       mon0 = mod(number_month-1,12) + 1
       imd = nmonth(mon0)
           if (mytid==0) write(6,*)'test-mm,imd,mon0,month',imd,mon0,month
       DO iday =number_day,imd 
       t3a=clock0f()
       !call COREIDEAL_DAILY
       if(1) call jra_daily
       t1=clock0f()
        if(mytid==0) write(*,*) 'jra_daily time :',t1-t3a
#ifndef GPU
       DO II = 1, NSS
#else GPU
       DO II = 1, NSS, NSS
#endif GPU

#endif COUP

#ifndef GPU
              if (II==1)call energy
         num_step_per_day=num_step_per_day+1
           curr_ymd_licom=iyfm*10000+mon0*100+iday   
!           if (mytid==0) write(6,*)'curr_ymd_cpl,tod_cpl',curr_ymd_cpl,tod_cpl
!---------------------------------------------------------------------
! set the ideal temperature and salinity for ideal simulation !LPF 20200304
!if the true value in the inirun, then this set will be undone!!
#if (defined TEST_IDEAL)
    if(0)  CALL setidealvalue
#endif

!---------------------------------------------------------------------
!     COMPUTE DENSITY, BAROCLINIC PRESSURE AND THE RELAVANT VARIABLES
    LOGMSG()
            CALL READYT

!---------------------------------------------------------------------
!     BAROCLINIC & BAROTROPIC CYCLE
!---------------------------------------------------------------------
    LOGMSG()
            DO JJ = 1,NCC

!     COMPUTE MOMENTUM ADVECTION, DIFFUSION & THEIR VERTICAL INTEGRALS
    LOGMSG()
            if (mytid ==0) then
                write(111,*)"OK----8.0"
                close(111)
            end if
               CALL READYC
            if (mytid ==0) then
                write(111,*)"OK----8.1"
                close(111)
            end if
!     PREDICTION OF BAROTROPIC MODE
    LOGMSG()
               CALL BAROTR
            if (mytid ==0) then
                write(111,*)"OK----8.2"
                close(111)
            end if
!     PREDICTION OF BAROCLINIC MODE
    LOGMSG()
               CALL BCLINC
            END DO

!--------------------------------------------------------------------
!     PREDICTION OF PASSIVE TRACER
!--------------------------------------------------------------------
    LOGMSG()
            CALL TRACER
            if (mytid ==0) then
                write(111,*)"OK----8.3"
                close(111)
            end if
     !   if(mytid==0)then
     !     write(*,*)"at1,tracer",at(2,2,26,1),at(2,2,27,1)
     !     write(*,*)"at2,tracer",at (2,2,26,2),at(2,2,27,2)
     !  endif
    LOGMSG()
            CALL ICESNOW
              call energy

      !  if(mytid==0)then
      !    write(*,*)"at1,ice",at(2,2,26,1),at(2,2,27,1)
      !    write(*,*)"at2,ice",at (2,2,26,2),at(2,2,27,2)
      ! endif

!--------------------------------------------------------------------
!     PERFORM CONVECTIVE ADJUSTMENT IF UNSTABLE STRATIFICATION OCCURS
!--------------------------------------------------------------------
    LOGMSG()
!       IF (mod(IST,(1800/IDTS))==0 ) CALL CONVADJ
        CALL CONVADJ
!--------------------------------------------------------------------
!Accumulate the data in a T/S time step interval or in a day for output
          ! if (mytid==0) write(6,*)'begin accumm'
             if (dts_accum==.true.) CALL ACCUMM
!--------------------------------------------------------------------
!       Finish thermal cycle loop in a couple interval
       if(mod(II,72) ==0) then
        t2=clock0f()
        if(mytid==0) write(*,*) ' (CPU) stepon test :',t2-t1, 'loops :',II,'/',NSS
        t1=clock0f()
       end if
#else
            call stepon(nss,1,ncc,0)
         num_step_per_day=NSS
        t2=clock0f()
        if(mytid==0) write(*,*) '(GPU/HIP) stepon test :',t2-t1, 'loops :',NSS
        t1=clock0f()
              call energy
#endif
   END DO 
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!          Diagnose if the energy is out of control
!--------------------------------------------------------------------

      if(num_step_per_day==NSS) then 
           num_step_per_day=0
           CALL ADDPS
      end if
!--------------------------------------------------------------------
!Accumulate the data in a T/S time step interval or in a day for output
          ! if (mytid==0) write(6,*)'begin accumm'
            if ((.not.dts_accum).and.(num_step_per_day==0)) CALL ACCUMM
           !if (mytid==0) write(6,*)'end accumm'
!--------------------------------------------------------------------

!     ACCUMULATE SOME VARIABLES FOR MONTHLY OUTPUT
         IF (iday ==imd .and.num_step_per_day==0) then
            month = month + 1
            dd_licom=1 !LPF20200315
            if (imt_global < 3000) CALL SSAVEMON !turnon onlyfor low
         ENDIF

#ifdef COUP
         call shr_sys_flush(6)
#endif

!--------------------------------------------------------------------
!     SAVE FOR OUTPUT DAILY
    LOGMSG()
     if ( num_step_per_day== 0) then
     t2=clock0f()
         CALL SSAVEINS
     t3b=clock0f()
     if(mytid==0) write(*,*) 'write_file SSH time :',t3b-t2
         CALL NEXTSTEP
     end if
        !if (mytid==0) write(6,*)'end SSAVE-INST'
!--------------------------------------------------------------------
!#if !defined(COUP)
       t3b=clock0f()
       if(mytid==0) write(*,FMT='(A6,F20.15,A10,F20.15,A6, I2,I3)') ' day:',t3b-t3a,' SYPD:',86400./((t3b-t3a)*365.),'   ',MM,iday
    end do
    end do
!#endif 
!--------------------------------------------------------------------
     call mpi_finalize(ierr)
     end 
