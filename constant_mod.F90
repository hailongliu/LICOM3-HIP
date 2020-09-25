
  module constant_mod

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use pmix_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
use shr_sys_mod
!use shr_msg_mod 
use shr_kind_mod 
use shr_const_mod
use netcdf

  implicit none

  public  :: const
!
  !----- constants -----
  real(SHR_KIND_R8), parameter ::  latvap   = SHR_CONST_LATVAP ! latent heat of evap   ~ J/kg
  real(SHR_KIND_R8), parameter ::  Tzro     = SHR_CONST_TKFRZ  ! 0 degrees C                       ~ kelvin
  real(SHR_KIND_R8), parameter ::  Tfrz     = Tzro   - 1.8     ! temp of saltwater freezing ~ kelvin
! real(SHR_KIND_R8), parameter ::  pi       = SHR_CONST_PI     ! a famous math constant
  real(SHR_KIND_R8), parameter ::  pi       = 4.0*ATAN(1.0)
! real(SHR_KIND_R8), parameter ::  omega    = SHR_CONST_OMEGA  ! earth's rotation  ~ rad/sec
  real(SHR_KIND_R8), parameter ::  omega    = 0.7292D-4
! real(SHR_KIND_R8), parameter ::  g        = SHR_CONST_G      ! gravity ~ m/s^2
  real(SHR_KIND_R8), parameter ::  g        = 9.806D0
  real(SHR_KIND_R8), parameter ::  DEGtoRAD = PI/180.0         ! PI/180
  real(SHR_KIND_R8), parameter ::  RADIUS   = 6371000D0
  real(r8), parameter :: very_small   =   1.d-15
  real(r8), parameter :: karman = 0.4d0


!
   real (r8), parameter, public :: &
      c0     =    0.0_r8   ,&
      c1     =    1.0_r8   ,&
      c2     =    2.0_r8   ,&
      c3     =    3.0_r8   ,&
      c4     =    4.0_r8   ,&
      c5     =    5.0_r8   ,&
      c8     =    8.0_r8   ,&
      c10    =   10.0_r8   ,&
      c16    =   16.0_r8   ,&
      c1000  = 1000.0_r8   ,&
      c10000 =10000.0_r8   ,&
      c1p5   =    1.5_r8   ,&
      p33    = c1/c3       ,&
      p5     = 0.500_r8    ,&
      p25    = 0.250_r8    ,&
      p125   = 0.125_r8    ,&
      p001   = 0.001_r8    ,&
      eps    = 1.0e-10_r8  ,&
      eps2   = 1.0e-20_r8  ,&
      bignum = 1.0e+30_r8  ,&
      pi2 = c2*pi


   real (r4), parameter, public ::         &
      undefined_nf_r4  = NF90_FILL_FLOAT,  &
      undefined        = -12345._r4

   real (r8), parameter, public ::         &
      undefined_nf_r8  = NF90_FILL_DOUBLE

   real (r8), public ::  &
      undefined_nf = NF90_FILL_DOUBLE

   integer (int_kind), parameter, public ::   &
      undefined_nf_int = NF90_FILL_INT

   !*** location of fields for staggered grids

   integer (int_kind), parameter, public ::   &
      field_loc_unknown  =  0, &
      field_loc_noupdate = -1, &
      field_loc_center   =  1, &
      field_loc_SWcorner =  2, &
      field_loc_Sface    =  3, &
      field_loc_Wface    =  4

   !*** field type attribute - necessary for handling
   !*** changes of direction across tripole boundary

   integer (int_kind), parameter, public ::   &
      field_type_unknown  =  0, &
      field_type_noupdate = -1, &
      field_type_scalar   =  1, &
      field_type_vector   =  2, &
      field_type_angle    =  3
!
  character (5), parameter, public :: &
      blank_fmt = "(' ')"

!
      contains 
!  CVS: $Id: const.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ================
      SUBROUTINE CONST
!     ================
!-----------------------------------------------------------------------
!
! Purpose: Set up some constants and control parameter.
!
! Author: Yongqiang Yu and Hailong Liu, Dec. 31, 2002
!
!-----------------------------------------------------------------------

!lhl090729
      integer :: ncid,iret
      real(r8) :: tmpy(s_jmt)
!lhl090729
      INTEGER :: IDTB,IDTC,IDTS
      REAL(r8)    :: AG,ALFA

      namelist /namctl/ AFB1,AFC1,AFT1,IDTB,IDTC,IDTS,AMV,AHV,NUMBER, &
                        NSTART,IO_HIST,IO_REST,klv,kvt,AM_HOR,AH_HOR,   &
                        hist_freq,out_dir,adv_tracer,adv_momentum, &
                        diag_msf,diag_bsf,diag_budget,test_input,diag_mth,rest_freq, boundary_restore, &
                        AM_BIHAR, AH_BIHAR,GAMMA,num_cpl,simple_assm,dailybudget_accum,&
                        dts_accum,daily_accum,yearadd,dayadd,iyfmfnew
!-------------------------------------------------------
!     Set up the default value for namelist.
!-------------------------------------------------------
!
!      diag_msf=.false.
!      diag_bsf=.false.
!      diag_budget=.false.
!      diag_mth=.false.

!

      dts_accum = .false.
      !dts_accum = .true.!LPF20200321
      daily_accum=.true.
      dailybudget_accum=.true.
      simple_assm=.false.


      IDTB=30   ; IDTC=1800 ; IDTS=7200
      AFB1=0.025D0; AFC1=0.43D0 ; AFT1=0.43D0
      num_cpl=1
      iyfmfnew=0
!-------------------------------------------------------
!     Set up the constants for calendar month.
!-------------------------------------------------------
      NMONTH=RESHAPE((/31,28,31,30,31,30,31,31,30,31,30,31/),(/12/))
     NNMONTH=RESHAPE((/0,31,59,90,120,151,181,212,243,273,304,334/),(/12/))
      ABMON =RESHAPE((/'Jan','Feb','Mar','Apr','May','Jun', &
                   'Jul','Aug','Sep','Oct','Nov','Dec'/),(/12/))
      ABMON1=RESHAPE((/'jan','feb','mar','apr','may','jun', &
                    'jul','aug','sep','oct','nov','dec'/),(/12/))
!-------------------------------------------------------
!     PHYSICAL CONSTANTS
!-------------------------------------------------------
!     CP    Specific heat capacity of sea water in J/kg/K
!     D0    Density of sea water in kg/m**3
!     AG    Ekman bias angle
!     C0F    friction coefficient
!     TBICE the frozen point of seawater in C

      CP = 3996.0D0
      D0 = 1026.0D0
      AG = 10.0D0*3.1415926D0/180.0D0
      C0F = 2.6D-3
      TBICE = -1.8D0

      SAG = SIN (AG)
      CAG = COS (AG)
      OD0 = 1.0D0/ D0
      OD0CP = 1.0D0/ (D0* CP)
      kvt=45 !the upper 1500m if total Layer 55
!-------------------------------------------------------
!     DIFFUSION & VISCOSITY
!-------------------------------------------------------
!     AM    laternal viscosity coeffcient in m**2/sec
!     AH    laternal diffusion coeffcient in m**2/sec
!     AMV   vertical viscosity coeffcient in m**2/sec
!     AHV   vertical diffusion coeffcient in m**2/sec
!     AHICE diffusion coeffcient between ice & water

      AMV = 1.0D-3
!lhl  AH  = 2.0D+3
      AHV = 0.3D-4

!
!-------------------------------------------------------
!     Read namelist from the external file.
!-------------------------------------------------------
      if (mytid==0) then
	      open(11,file='ocn.parm',form='formatted',status='OLD')
	      rewind 11
	      read(11,namctl)
	      close(11)
      end if
!     call mpi_barrier(mpi_comm_ocn,ierr)
      call mpi_bcast(afb1,1,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(afc1,1,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(aft1,1,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(amv,1,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(ahv,1,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(idtb,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(idtc,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(idts,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(number,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(nstart,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(io_hist,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(io_rest,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(klv,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(kvt,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(am_hor,1,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(ah_hor,1,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(am_bihar,1,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(ah_bihar,1,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(hist_freq,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(rest_freq,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(adv_tracer,80,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(adv_momentum,80,mpi_character,0,mpi_comm_ocn,ierr)
      call mpi_bcast(out_dir,80,mpi_character,0,mpi_comm_ocn,ierr)
      call mpi_bcast(nstart,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(boundary_restore,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(diag_msf,1,mpi_logical,0,mpi_comm_ocn,ierr)
      call mpi_bcast(diag_bsf,1,mpi_logical,0,mpi_comm_ocn,ierr)
      call mpi_bcast(diag_budget,1,mpi_logical,0,mpi_comm_ocn,ierr)
      call mpi_bcast(daily_accum,1,mpi_logical,0,mpi_comm_ocn,ierr)
      call mpi_bcast(dailybudget_accum,1,mpi_logical,0,mpi_comm_ocn,ierr)
      call mpi_bcast(dts_accum,1,mpi_logical,0,mpi_comm_ocn,ierr)
      call mpi_bcast(simple_assm,1,mpi_logical,0,mpi_comm_ocn,ierr)
      call mpi_bcast(test_input,1,mpi_logical,0,mpi_comm_ocn,ierr)
      call mpi_bcast(diag_mth,1,mpi_logical,0,mpi_comm_ocn,ierr)
      call mpi_bcast(gamma,1,mpi_pr,0,mpi_comm_ocn,ierr)
      call mpi_bcast(num_cpl,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(yearadd,1,mpi_integer,0,mpi_comm_ocn,ierr) !LPF20170621
      call mpi_bcast(dayadd,1,mpi_integer,0,mpi_comm_ocn,ierr)  !LPF20170621
      call mpi_bcast(iyfmfnew,1,mpi_integer,0,mpi_comm_ocn,ierr)  !LPF20170621
!      call mpi_barrier(mpi_comm_ocn,ierr)

      AHICE = AHV

!-------------------------------------------------------
!     wndmix = min value for mixing at surface to simulate
!              high freq wind mixing. (m**2/sec)
!     fricmx = maximum mixing (m**2/sec)
!     diff_cbt_back = background "diff_cbt"(m**2/sec)
!     visc_cbu_back = background"visc_cbu" (m**2/sec)
!     diff_cbt_limit = largest "diff_cbt"  (m**2/sec)
!     visc_cbu_limit = largest "visc_cbu"  (m**2/sec)
!-------------------------------------------------------

      wndmix = 10.0d-4
      fricmx = 50.0d-4
      diff_cbt_back = 1.0d-4
      visc_cbu_back = 1.0d-4

      visc_cbu_limit = fricmx
      diff_cbt_limit = fricmx

!lhl1204
! 1 in 2x2 model, baroclinic time step is 7200s at first and 9600s later,
!   and annual mean forcina, and from the begining
!      dfricmx = 1000.0d-4  !OK
!   more larger parameter is tested!
!   the time step of baroclinic and thermal have been reduced
!   to 1/4 (3600) and 1/2 (14400) respectively.
!   it is same as when the viscosity decreased.
!      dfricmx = 2000.0d-4  !OK
! 2 this parameter can also use when the viscosity decreased to 1e+4
!   but the time step of baroclinic and thermal have been reduced to
!   1/4 (3600) and 1/2 (14400) respectively. the former is close to the
!   time step which used in Canoto's 2001 paper.
!      dfricmx = 1000.0d-4  !OK
! in 2x2 model, baroclinic time step is 7200s, seasonal forcing, from the first experiment
!      dfricmx = 500.0d-4  !OK
       dfricmx = 100.0d-4 !OK
!       dfricmx = 2000.0d-4 !OK
!      dfricmx = 5000.0d-4 !OK
      dwndmix =   10.0d-4
!lhl1204

!LPF20160505
#ifdef BCKMEX
   diff_back_coef_max = 0.13
   diff_back_eq       = 0.01
   diff_back_coef     = 0.16
#endif
   Pr_number=10.0
!LPF20160505


      DTB = FLOAT (IDTB)
      DTC = FLOAT (IDTC)
      DTS = FLOAT (IDTS)

      DTB2 = 2.0D0* DTB
      DTC2 = 2.0D0* DTC

!-------------------------------------------------------
!     NSS   steps for thermohaline procoss per day
!     NCC   baroclinic steps within one thermohaline step
!     NBB   barotropic steps within one baroclinic step
!-------------------------------------------------------
      NSS = 24*3600/ IDTS
#if ( defined SYNCH)
      NCC = IDTS / IDTC
#else
      NCC = 1
#endif
      NBB = IDTC / IDTB

      ONBB = 1.0D0/ FLOAT (NBB +1)
      ONCC = 1.0D0/ FLOAT (NCC +1)
      ONBC = 1.0D0/ FLOAT (NBB * NCC +1)

      AFB2 = 1.0D0-2.0D0* AFB1
      AFC2 = 1.0D0-2.0D0* AFC1
      AFT2 = 1.0D0-2.0D0* AFT1


!-------------------------------------------------------
!     RESTORING TIME SCALE FOR SALINITY
!-------------------------------------------------------
      GAMMA = 1.0D0/ (GAMMA*86400.0D0)
!YU

!-------------------------------------------------------
!     Read reference temperature, salinity and
!     coefficients for UNESCO formula.
!-------------------------------------------------------
      if (mytid==0) then
      open(33,file='dncoef.h1',form='formatted')
      read(33,*)
      read(33,*)to
      read(33,*)
      read(33,*)so
!lhl1204   add a arry to read reference potential density
!
      read(33,*)
      read(33,*)po
!      read(33,'(8(f7.4,1x))')po
!lhl1204
      do i=1,km
         read(33,*)
         read(33,*)(c(i,k),k=1,9)
      end do
      close(33)
      endif
!     call mpi_barrier(mpi_comm_ocn,ierr)
      call mpi_bcast(to,km,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(so,km,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(po,km,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(c,km*9,MPI_PR,0,mpi_comm_ocn,ierr)
!     call mpi_barrier(mpi_comm_ocn,ierr)
!YU
!-------------------------------------------------------
!     the units of density departures is [g/cm**3] in
!     GFDL's MOM, and [kg/m**3] in present model
!-------------------------------------------------------

!$OMP PARALLEL DO PRIVATE (M,K)
      DO M = 1,9
         DO K = 1,KM
            C (K,M)= C (K,M)*1000.0D0
         END DO
      END DO

!YU
      if (mytid ==0 ) then
      write(6,*)"AM_HOR,AH_HOR",AM_HOR,AH_HOR
      endif
      !dts_accum = .false.
     ! dts_accum = .true.
     ! daily_accum=.true.
     ! dailybudget_accum=.true.
     ! simple_assm=.true.
!YU

      RETURN
      END SUBROUTINE CONST
!
     end module constant_mod
