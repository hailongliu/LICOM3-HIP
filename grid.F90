!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module grid

!BOP
! !MODULE: grid
!
! !DESCRIPTION:
!  This module contains grid info and routines for setting up the
!  POP grid quantities.
!
! !REVISION HISTORY:
!  SVN:$Id: grid.F90 17212 2009-07-20 23:01:42Z njn01 $

! !USES:

   use precision_mod
   use param_mod
   use LICOM_Error_mod
   use POP_HaloMod
   use POP_GridHorzMod
   use blocks
   use distribution
   use domain
   use broadcast
   use gather_scatter
   use constant_mod
   use msg_mod
   use global_reductions
   use pconst_mod
   use cdf_mod !for reading basin index filed

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public  :: init_grid1,     &
              init_grid2,     &
              read_horiz_grid,&
              tgrid_to_ugrid, &
              ugrid_to_tgrid, &
              read_topography,&
              calc_coeff

! !PUBLIC DATA MEMBERS:

   real (r8), public :: &
      area_u, area_t       ,&! total ocean area of U,T cells
      volume_u, volume_t   ,&! total ocean volume of U,T cells
      volume_t_marg        ,&! volume of marginal seas (T cells)
      area_t_marg          ,&! area of marginal seas (T cells)
      uarea_equator          ! area of equatorial cell

   real (r8), dimension(km), public :: &
      area_t_k             ,&! total ocean area (T cells) at each dpth
      volume_t_k           ,&! total ocean volume (T cells) at each dpth
      volume_t_marg_k        ! tot marginal seas vol (T cells) at each dpth

   integer (i4), public :: &
      sfc_layer_type,       &! choice for type of surface layer
      kmt_kmin, 	    &! minimum allowed non-zero KMT value
      n_topo_smooth 	     ! number of topo smoothing passes

   integer (i4), parameter, public :: &
      sfc_layer_varthick = 1,  &! variable thickness surface layer
      sfc_layer_rigid    = 2,  &! rigid lid surface layer
      sfc_layer_oldfree  = 3    ! old free surface form

   logical (log_kind), public ::    &
      partial_bottom_cells   ! flag for partial bottom cells

   real (r8), dimension(:,:), allocatable, public :: &
      BATH_G           ! Observed ocean bathymetry mapped to global T grid
                       ! for use in computing KMT internally

   integer (i4), dimension(:,:), allocatable, public :: &
      KMT_G            ,&! k index of deepest grid cell on global T grid
      BASIN_G   !zwp         ! for use in performing work distribution


!-----------------------------------------------------------------------
!
!  grid information for all local blocks
!  the local blocks are by default in baroclinic distribution
!
!-----------------------------------------------------------------------


   !*** geometric 2d arrays

   real (r8), allocatable,dimension(:,:,:), public :: &
      DXU, DYU            ,&! {x,y} spacing centered at U points
      DXT, DYT            ,&! {x,y} spacing centered at T points
      DXUR, DYUR          ,&! reciprocals of DXU, DYU
      DXTR, DYTR          ,&! reciprocals of DXT, DYT
      HTS, HTW            ,&! cell widths on {S,W} sides of T cell
      HUN, HUE            ,&! cell widths on {N,E} sides of U cell
      ULAT, ULON          ,&! {latitude,longitude} of U points
      TLAT, TLON          ,&! {latitude,longitude} of T points
      ANGLE, ANGLET       ,&! angle grid makes with latitude line
      FCOR, FCORT         ,&! coriolis parameter at U,T points
      UAREA, TAREA        ,&! area of U,T cells
      UAREA_R, TAREA_R    ,&! reciprocal of area of U,T cells
      HT, HU, HUR           ! ocean depth at T,U points

   !*** 3d depth fields for partial bottom cells

   real (r8), dimension(:,:,:,:), allocatable, public :: &
      DZU, DZT               ! thickness of U,T cell for pbc

   !*** 2d landmasks
!ZWP

   integer, dimension(imt,jmt,max_blocks_clinic), &
      public :: &
      BASIN            ! basin index of deepest grid cell on T grid
!ZWP

   integer (i4), allocatable,dimension(:,:,:), &
      public :: &
      KMT            ,&! k index of deepest grid cell on T grid
      KMU            ,&! k index of deepest grid cell on U grid
      KMTOLD           ! KMT field before smoothing

   logical (log_kind), dimension(nx_block,ny_block,max_blocks_clinic), &
      public :: &
      CALCT          ,&! flag=true if point is an ocean point
      CALCU            !   at the surface

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), public :: &
      RCALCT         ,&! real equiv of CALCT,U to use as more
      RCALCU           !   efficient multiplicative mask

   integer (i4), allocatable,dimension(:,:,:), &
      public :: &
      KMTN,KMTS,KMTE,KMTW   ,&! KMT field at neighbor points
      KMUN,KMUS,KMUE,KMUW     ! KMU field at neighbor points

   integer (i4), dimension(nx_block,ny_block,max_blocks_clinic), &
      public :: &
      KMTEE,KMTNN      ! KMT field 2 cells away for upwind stencil
                       ! allocated and computed in advection module

   integer (i4), dimension(:,:,:), allocatable, public :: &
      REGION_MASK      ! mask defining regions, marginal seas

!-----------------------------------------------------------------------
!
!     define types used with region masks and marginal seas balancing
!
!-----------------------------------------------------------------------
   integer (i4), parameter, public :: &
      max_regions =   15, &              ! max no. ocean regions
      max_ms      =    7                 ! max no. marginal seas
 
   integer (i4), public :: &
      num_regions, &
      num_ms

   type, public :: ms_bal
      real    (r8)         :: lat         ! transport latitude
      real    (r8)         :: lon         ! transport longitude
      real    (r8)         :: area        ! total distribution area
      real    (r8)         :: transport   ! total excess/deficit (E+P+M+R)
      integer (i4)   :: mask_index  ! index of m-s balancing mask
   end type ms_bal

   type, public :: regions                ! region-mask info
      integer   (i4) :: number 
      character (char_len) :: name          
      logical   (log_kind) :: marginal_sea
      real      (r8      ) :: area
      real      (r8      ) :: volume
      type      (ms_bal)   :: ms_bal
   end type regions       

   type (regions),dimension(max_regions), public :: region_info

   integer (i4),public ::       &
      nocean_u, nocean_t,      &! num of ocean U,T points
      nsurface_u, nsurface_t    ! num of ocean U,T points at surface

!#################### temporary kludge for overflows ####################
   character (char_len), public ::  &
      topography_filename    ! copy of input topography filename 
!########################################################################

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  module private data
!
!-----------------------------------------------------------------------

   !*** geometric scalars

   integer (i4) ::       &
      jeq                       ! j index of equatorial cell

   logical (log_kind) ::  &
      flat_bottom,        &! flag for flat-bottom topography
      lremove_points       ! flag for removing isolated points

   real (r8), dimension(:,:), allocatable :: &
      TLAT_G, TLON_G        ! {latitude,longitude} of U points
                            ! in global-sized array

!LPF20160819
   real (r8), dimension(:,:), allocatable :: &
      ULAT_G, ULON_G        ! {latitude,longitude} of U points
!LPF20160819
!-----------------------------------------------------------------------
!
!     area-weighted averaging coefficients
!     AT{0,S,W,SW} = {central,s,w,sw} coefficients for area-weighted
!       averaging of four U points surrounding a T point
!     AU{0,N,E,NE} = {central,n,e,ne} coefficients for area-weighted
!       averaging of four T points surrounding a U point
!
!-----------------------------------------------------------------------

   real (r8), allocatable,dimension (:,:,:), public :: &
      AT0,ATN,ATE,ATNE,AU0,AUS,AUW,AUSW

!-----------------------------------------------------------------------
!
!  variables which are shared between init_grid1,init_grid2
!
!-----------------------------------------------------------------------

   character (char_len), public ::  &
      horiz_grid_opt,       &! horizontal grid option
      vert_grid_opt,        &! vertical grid option
      sfc_layer_opt,        &! choice for surface layer type
      topography_opt,       &! topography (KMT) option
      horiz_grid_file,      &! input file for reading horiz grid info
      vert_grid_file,       &! input file for reading horiz grid info
      topography_file,      &! input file for reading horiz grid info
      region_mask_file,     &! input file for region mask
      region_info_file,     &! input file with region identification
      bottom_cell_file,     &! input file for thickness of pbc
      topography_outfile,   &! output file for writing horiz grid info
      basin_grid_file        !  input file for basin index

    integer :: stdout
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_grid1
! !INTERFACE:

 subroutine init_grid1

! !DESCRIPTION:
!  Initializes only grid quantities necessary for completing
!  decomposition setup (ULAT, ULON, KMT).
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   namelist /grid_nml/ horiz_grid_file, vert_grid_file, topography_file, horiz_grid_opt, basin_grid_file

   integer (i4) :: &
      nml_error           ! namelist i/o error flag

!-----------------------------------------------------------------------
!
!  read input namelist for grid setup options
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  get global ULAT,ULON
!
!-----------------------------------------------------------------------
   if (my_task == master_task) then
      open (11, file='ocn.parm', status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(11, nml=grid_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(11)
      write(6,*) " horiz_grid_file, vert_grid_file, topography_file, horiz_grid_opt, basin_grid_file"
      write(6,*) horiz_grid_file, vert_grid_file, topography_file, horiz_grid_opt, basin_grid_file
      call flush(6)
   endif


      call broadcast_scalar(horiz_grid_file, master_task)
      call broadcast_scalar(vert_grid_file, master_task)
      call broadcast_scalar(topography_file, master_task)
      call broadcast_scalar(horiz_grid_opt, master_task)
      call read_horiz_grid(horiz_grid_file,.true.)

!-----------------------------------------------------------------------
!  set up topography by getting global KMT field (used for
!  creating a load balanced block distribution).
!
!-----------------------------------------------------------------------

      call broadcast_scalar(topography_file, master_task)
      call read_topography(topography_file,.true.)

 end subroutine init_grid1

!***********************************************************************
!BOP
! !IROUTINE: init_grid2
! !INTERFACE:

 subroutine init_grid2

! !DESCRIPTION:
!  Initializes all grid quantities
!
! !REVISION HISTORY:
!  same as module
!
! !OUTPUT PARAMETERS:

   integer (i4) :: errorCode              ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: &
      reclength,ioerr,nu,i,j,k,n,iblock,    &! dummy loop index variables
      range_count         ! counter for angle out of range

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      WORK                 ! local temp space

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      DZBC                 ! thickness of bottom T cell for pbc

   character (*), parameter ::  &! output formats
      vgrid_fmt1 = "(3x,' k ',3x,'Thickness (cm)',3x,' Depth (cm) ')", &
      vgrid_fmt2 = "(3x,'---',3x,'--------------',3x,'------------')", &
      vgrid_fmt3 = "(3x,i3,4x,1pe12.5,4x,1pe12.5)"                   , &
      topo_fmt1  = "(' # surface (T,U) points',2x,i10,2x,i10)"       , &
      topo_fmt2  = "(' # ocean   (T,U) points',2x,i10,2x,i10)"       , &
      topo_fmt3  = "(' T-area,   U-area   (km^2)',2(2x,1pe23.15))"   , &
      topo_fmt4  = "(' T-volume, U-volume (km^3)',2(2x,1pe23.15))"

   type (block) :: &
      this_block  ! block info for current block

   real (r8) ::         &
      angle_0, angle_e, &! temporaries for computing angle at T points
      angle_n, angle_ne 

!-----------------------------------------------------------------------
!
!  output grid setup options to log file
!
!-----------------------------------------------------------------------

   errorCode = 0

   stdout = 6
   if (my_task == master_task) then
      write(stdout,'(a13)') ' Grid options'
   endif

!-----------------------------------------------------------------------
!
!  set up horizontal grid
!
!-----------------------------------------------------------------------

      call broadcast_scalar(horiz_grid_file, master_task)
      call read_horiz_grid(horiz_grid_file,.false.)


!-----------------------------------------------------------------------
!
!  if boundaries are closed, extend physical domain values into ghost
!  cells
!  compute other derived quantities like areas and reciprocals
!
!-----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE(iblock,i,j,this_block)
   do iblock=1,nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)

      if (this_block%i_glob(1) == 0) then ! closed western bndy
         do j=1,ny_block
         do i=1,this_block%ib-1
            DXU(i,j,iblock) = DXU(this_block%ib,j,iblock)
            DYU(i,j,iblock) = DYU(this_block%ib,j,iblock)
            DXT(i,j,iblock) = DXT(this_block%ib,j,iblock)
            DYT(i,j,iblock) = DYT(this_block%ib,j,iblock)
         end do
         end do
      endif

      if (this_block%i_glob(this_block%ie+1) == 0) then ! closed east bndy
         do j=1,ny_block
         do i=this_block%ie+1,nx_block
            DXU(i,j,iblock) = DXU(this_block%ie,j,iblock)
            DYU(i,j,iblock) = DYU(this_block%ie,j,iblock)
            DXT(i,j,iblock) = DXT(this_block%ie,j,iblock)
            DYT(i,j,iblock) = DYT(this_block%ie,j,iblock)
         end do
         end do
      endif

      if (this_block%j_glob(1) == 0) then ! closed north bndy
         do j=1,this_block%jb-1
         do i=1,nx_block
            DXU(i,j,iblock) = DXU(i,this_block%jb,iblock)
            DYU(i,j,iblock) = DYU(i,this_block%jb,iblock)
            DXT(i,j,iblock) = DXT(i,this_block%jb,iblock)
            DYT(i,j,iblock) = DYT(i,this_block%jb,iblock)
         end do
         end do
      endif

      if (this_block%j_glob(this_block%je+1) == 0) then ! closed south bndy
         do j=this_block%je+1,ny_block
         do i=1,nx_block
            DXU(i,j,iblock) = DXU(i,this_block%je,iblock)
            DYU(i,j,iblock) = DYU(i,this_block%je,iblock)
            DXT(i,j,iblock) = DXT(i,this_block%je,iblock)
            DYT(i,j,iblock) = DYT(i,this_block%je,iblock)
         end do
         end do
      endif

      DXUR(:,:,iblock) = c1/DXU(:,:,iblock)
      DYUR(:,:,iblock) = c1/DYU(:,:,iblock)

      UAREA(:,:,iblock) = DXU(:,:,iblock)*DYU(:,:,iblock)
      UAREA_R(:,:,iblock) = c1/UAREA(:,:,iblock)

      DXTR(:,:,iblock) = c1/DXT(:,:,iblock)
      DYTR(:,:,iblock) = c1/DYT(:,:,iblock)

      TAREA(:,:,iblock) = DXT(:,:,iblock)*DYT(:,:,iblock)
      TAREA_R(:,:,iblock) = c1/TAREA(:,:,iblock)

   end do
!!$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!  calculate stencil coefficients for area-averaging
!
!-----------------------------------------------------------------------

   call cf_area_avg  ! coefficients for area-weighted averages

!-----------------------------------------------------------------------
!
!  calculate lat/lon of T points and calculate ANGLET from ANGLE
!
!-----------------------------------------------------------------------

   call calc_upoints(errorCode)

   if (errorCode /= 0) then
      call LICOM_ErrorSet(errorCode, &
         'init_grid2: error in calc_upoints')
      return
   endif

!LPF20160819
!for output the lon/lat on u-grid
      if (.not. allocated(ULAT_G)) then
         allocate (ULAT_G(imt_global,jmt_global), &
                   ULON_G(imt_global,jmt_global))
      endif
   call  gather_global(ULAT_G, ULAT, master_task,distrb_clinic) !LPF20160819
   call  gather_global(ULON_G, ULON, master_task,distrb_clinic) !LPF20160819

    ulat_o=ULAT_G/DEGtoRAD
    ulon_o=ULON_G/DEGtoRAD

    deallocate(ULAT_G,ULON_G)
!LPF20160819

!  !***
!  !*** first, ensure that -pi <= ANGLET <= pi
!  !***

   range_count = global_count ((ANGLET < - pi .or. ANGLET > pi), &
                               distrb_clinic, field_loc_SWcorner)

   if (range_count > 0) call exit_LICOM(sigAbort, &
                        'ERROR: ANGLET is outside its expected range')

!  !***
!  !*** compute ANGLE on T-grid
!  !***

!   ANGLET = c0


   !$OMP PARALLEL DO PRIVATE (n,i,j,angle_0,angle_e,angle_n,angle_ne, &
   !$OMP                      this_block)

   do n=1,nblocks_clinic
      this_block = get_block(blocks_clinic(n),n)

      do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie

            angle_0  = ANGLET(i,  j  ,n)
            angle_e  = ANGLET(i+1,j  ,n)
            angle_n  = ANGLET(i,  j-1,n)
            angle_ne = ANGLET(i+1,j-1,n)

            if ( angle_0 < c0 ) then
               if ( abs(angle_e -angle_0) > pi ) &
                  angle_e  = angle_e  - pi2 
               if ( abs(angle_n -angle_0) > pi ) &
                  angle_n  = angle_n  - pi2 
               if ( abs(angle_ne-angle_0) > pi ) &
                  angle_ne = angle_ne - pi2 
            endif

            ANGLE(i,j,n) =  angle_0 *AU0 (i,j,n) + &
                             angle_e *AUW (i,j,n) + &
                             angle_n *AUS (i,j,n) + &
                             angle_ne*AUSW(i,j,n)

         enddo

         !***
         !*** set ANGLET to zero for all of (global) j=1 row 
         !*** (bottom row of ANGLET is not used, but is written to file)
         !***

         if (this_block%j_glob(j) == 1) ANGLET(:,j,n) = c0

      enddo
   enddo
   !$OMP END PARALLEL DO

   call POP_HaloUpdate(ANGLE, POP_haloClinic, POP_gridHorzLocCenter, &
                               POP_fieldKindAngle, errorCode,         &
                               fillValue = 0.0_r8)

   if (errorCode /= 0) then
      call LICOM_ErrorSet(errorCode, &
         'init_grid2: error updating angleT halo')
      return
   endif

!-----------------------------------------------------------------------
!
!  set up vertical grid
!
!-----------------------------------------------------------------------

      if (my_task == master_task) then
         write(stdout,*) ' Reading vertical grid from file:', &
                         trim(vert_grid_file)
      endif
      call broadcast_scalar(vert_grid_file, master_task)
      call read_vert_grid(vert_grid_file)

!  !***
!  !*** calculate other vertical grid quantities
!  !***
!
!-----------------------------------------------------------------------
!
!  set up basin
!
!-----------------------------------------------------------------------
      if ( diag_mth .or. diag_msf) then
        if (my_task == master_task) then
           write(stdout,*) ' Reading Basin Index from file:', trim(basin_grid_file)
        endif
        call broadcast_scalar(basin_grid_file, master_task)
        call read_basin(basin_grid_file)
      end if
!      write(899,*) basin
!      stop

!-----------------------------------------------------------------------
!
!  set up topography
!
!-----------------------------------------------------------------------

      if (my_task == master_task) write(stdout,'(a30,a)') &
         ' Reading topography from file:', trim(topography_file)
      call broadcast_scalar(topography_file, master_task)
      call read_topography(topography_file,.false.)

!-----------------------------------------------------------------------
!
!  landmasks
!
!-----------------------------------------------------------------------
      do n=1,nblocks_clinic
         HT (:,:,n) = c0
         HU (:,:,n) = c0
         HUR(:,:,n) = c0

         do k=1,km
            do j=1,ny_block
            do i=1,nx_block
               if (k == KMT(i,j,n)) HT(i,j,n) = -zkp(k+1)
               if (k == KMU(i,j,n)) then
                  HU (i,j,n) = -zkp(k+1)
                  HUR(i,j,n) = -c1/zkp(k+1)
               endif
            enddo
            enddo
         enddo
      enddo


   call landmasks

!-----------------------------------------------------------------------
!
!  calculate area, volume, # surface points, # ocean points
!
!-----------------------------------------------------------------------

   area_t   = global_sum(TAREA, distrb_clinic, field_loc_center, RCALCT)
   area_u   = global_sum(UAREA, distrb_clinic, field_loc_SWcorner, RCALCU)
   WORK = TAREA*HT
   volume_t = global_sum(WORK, distrb_clinic, field_loc_center, RCALCT)
   WORK = UAREA*HU
   volume_u = global_sum(WORK, distrb_clinic, field_loc_SWcorner, RCALCU)
   area_t_k(1) = area_t
   volume_t_k(1) = global_sum(TAREA*dzp(1), distrb_clinic, &
                              field_loc_center, RCALCT)
   do k=2,km
      WORK = merge(TAREA, c0, k <= KMT)
      area_t_k(k) = global_sum(WORK, distrb_clinic, field_loc_center)
      WORK = merge(TAREA*dzp(k), c0, k <= KMT)
      volume_t_k(k) = global_sum(WORK, distrb_clinic, field_loc_center)
   end do

   nsurface_t = global_count(RCALCT, distrb_clinic, field_loc_center)
   nsurface_u = global_count(RCALCU, distrb_clinic, field_loc_SWcorner)
   nocean_t = global_sum(KMT, distrb_clinic, field_loc_center)
   nocean_u = global_sum(KMU, distrb_clinic, field_loc_SWcorner)

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,topo_fmt1) nsurface_t, nsurface_u
      write(stdout,topo_fmt2) nocean_t, nocean_u
      write(stdout,topo_fmt3) area_t*1.0e-6_r8, &
                              area_u*1.0e-6_r8
      write(stdout,topo_fmt4) volume_t*1.0e-9_r8, &
                              volume_u*1.0e-9_r8
      write(stdout,blank_fmt)
   endif

!
!  compute coriolis parameter 2*omega*sin(true_latitude)
!
!-----------------------------------------------------------------------

   FCOR  = c2*omega*sin(ULAT)    ! at u-points
   FCORT = c2*omega*sin(TLAT)    ! at t-points

!-----------------------------------------------------------------------
!EOC


 end subroutine init_grid2

!***********************************************************************
! !IROUTINE: read_horiz_grid
! !INTERFACE:

 subroutine read_horiz_grid(horiz_grid_file, latlon_only)

! !DESCRIPTION:
!  Reads horizontal grid information from input grid file
!
! !REVISION HISTORY:
!  same as module

!  !INPUT PARAMETERS:

   character (*), intent(in) :: &
      horiz_grid_file     ! filename of file containing grid data

   logical (log_kind), intent(in) :: &
      latlon_only       ! flag requesting only ULAT, ULON

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: &
      i,j,iblock        ,&! loop counters
      ip1, im1, jp1, jm1,&! shift indexes
      nu                ,&! i/o unit number
      ioerr             ,&! i/o error flag
      reclength           ! record length
   integer (i4) :: errorCode

   type (block) :: &
      this_block

   real (r8) :: temp_xxx(jmt_global)
!-----------------------------------------------------------------------
!
!  if only lat,lon are requested, read only these
!
!-----------------------------------------------------------------------

   if (latlon_only) then
   
      allocate (TLAT_G(imt_global,jmt_global), &
                TLON_G(imt_global,jmt_global)) 
                
      INQUIRE(iolength=reclength) TLAT_G
      if (my_task == master_task) then
         open(25,file=trim(horiz_grid_file),status='old', &
              form='unformatted', access='direct', recl=reclength, &
              iostat=ioerr)
      endif   
      
      call broadcast_scalar(ioerr, master_task)
      if (ioerr /= 0) call exit_licom(sigAbort, &
                                    'Error opening horiz_grid_file')
                                    
      if (my_task == master_task) then
         read(25,rec=1,iostat=ioerr) TLAT_G
         read(25,rec=2,iostat=ioerr) TLON_G
         close(25)
      endif
      
      call broadcast_scalar(ioerr, master_task)
      if (ioerr /= 0) call exit_LICOM(sigAbort, &
                                    'Error reading horiz_grid_file')
                                    
      call broadcast_array(TLAT_G, master_task)
      call broadcast_array(TLON_G, master_task)

!-----------------------------------------------------------------------
!
!  otherwise, read everything else
!  compute some derived fields here to preserve information that is
!    lost once land blocks are dropped
!
!-----------------------------------------------------------------------

   else

      if (.not. allocated(TLAT_G)) then
         allocate (TLAT_G(imt_global,jmt_global), &
                   TLON_G(imt_global,jmt_global))
      endif

      INQUIRE(iolength=reclength) TLAT_G

      if (my_task == master_task) then
         open(25,file=trim(horiz_grid_file),status='old', &
              form='unformatted', access='direct', recl=reclength,iostat=ioerr)
         read(25,rec=1,iostat=ioerr) TLAT_G
         read(25,rec=2,iostat=ioerr) TLON_G
      endif
!ZWP20131030
      lon_o = TLON_G/DEGtoRAD
      lat_o = TLAT_G/DEGtoRAD
!ZWP20131030
      do j=1, jmt_global
         lat(j) = tlat_g(1,j)/DEGtoRAD
      end do
      do i=1, imt_global
         lon(i) = tlon_g(i,1)/DEGtoRAD
      end do
!

      call scatter_global(TLAT, TLAT_G, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
      call scatter_global(TLON, TLON_G, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)

      if (my_task == master_task) then
         read(25,rec=3,iostat=ioerr) TLAT_G  ! holds HUN
      endif
!
     if(my_task == master_task .and. trim(horiz_grid_opt) == 'lat_lon') then
        open(35,file="DXT.DAT",form="unformatted")
        read(35) temp_xxx
        close(35)
      do j=1,jmt_global
         do i=1,imt_global
            TLAT_G(i,j) = temp_xxx(j)
         end do
      end do
     end if


      call scatter_global(HUN, TLAT_G, master_task, distrb_clinic, &
                          field_loc_Sface, field_type_scalar)

      do j=1,jmt_global
      do i=1,imt_global
         ip1 = i+1
         if (i == imt_global) ip1 = 1 ! assume cyclic. non-cyclic
                                     ! will be handled during scatter !DXU
         TLON_G(i,j) = p5*(TLAT_G(i,j) + TLAT_G(ip1,j))
      end do
      end do
      call scatter_global(DXT, TLON_G, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)

      do j=1,jmt_global
         jp1 = j+1
         if (j == jmt_global) jp1= jmt_global
         do i=1,imt_global
            !DXT = p5(HTS(i,j)+HTS(i,j-1))
            TLON_G(i,j) = p5*(TLAT_G(i,j) + TLAT_G(i,jp1))
         end do
      end do
!
     if(my_task == master_task .and. trim(horiz_grid_opt) == 'lat_lon') then
        open(35,file="DXU.DAT",form="unformatted")
        read(35) temp_xxx
        close(35)
      do j=1,jmt_global
         do i=1,imt_global
            TLON_G(i,j) = temp_xxx(j)
         end do
      end do
     end if
!

      call scatter_global(DXU, TLON_G, master_task, distrb_clinic, &
                          field_loc_swcorner, field_type_scalar)

      if (my_task == master_task) then
         read(25,rec=4,iostat=ioerr) TLAT_G  ! holds HUE
      endif

     if(my_task == master_task .and. trim(horiz_grid_opt) == 'lat_lon') then
        open(35,file="DYU.DAT",form="unformatted")
        read(35) temp_xxx
        close(35)
      do j=1,jmt_global
         do i=1,imt_global
            TLAT_G(i,j) = temp_xxx(j)
         end do
      end do
     end if
!
      call scatter_global(HUE, TLAT_G, master_task, distrb_clinic, &
                          field_loc_Wface, field_type_scalar)

      do j=1,jmt_global
      do i=1,imt_global
         im1 = i-1
         if (i == 1 ) im1 = imt_global ! assume cyclic. non-cyclic
                                     ! will be handled during scatter !DYT
         TLON_G(i,j) = p5*(TLAT_G(i,j) + TLAT_G(im1,j))
      end do
      end do
!
      call scatter_global(DYU, TLON_G, master_task, distrb_clinic, &
                          field_loc_swcorner, field_type_scalar)

      do j=1,jmt_global
         jm1 = j-1
         if (j == 1 ) jm1 = 1 ! assume cyclic. non-cyclic ! will be handled during scatter
         do i=1,imt_global

            !DYU = p5(HTW(i,j)+HTS(i,j+1))
            TLON_G(i,j) = p5*(TLAT_G(i,j) + TLAT_G(i,jm1))
         end do
      end do
!
!-----------------------------------------------------------------------
!   Add tripole-grid correction (KL and GD)
!
!-----------------------------------------------------------------------
      if (ltripole_grid) then
         j= 1
         do i=1,imt_global
            TLON_G(i,j) = TLAT_G(i,j)
         end do
      endif
!
     if(my_task == master_task .and. trim(horiz_grid_opt) == 'lat_lon') then
        open(35,file="DYT.DAT",form="unformatted")
        read(35) temp_xxx
        close(35)
      do j=1,jmt_global
         do i=1,imt_global
            TLON_G(i,j) = temp_xxx(j)
         end do
      end do
     end if

      call scatter_global(DYT, TLON_G, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)

      if (my_task == master_task) then
         read(25,rec=5,iostat=ioerr) TLAT_G
         read(25,rec=6,iostat=ioerr) TLON_G
      endif

      call scatter_global(HTS, TLAT_G, master_task, distrb_clinic, &
                          field_loc_Sface, field_type_scalar)
      call scatter_global(HTW, TLON_G, master_task, distrb_clinic, &
                          field_loc_Wface, field_type_scalar)
!
      if (trim(horiz_grid_opt) == 'lat_lon') then
         hun = dxt
         hue = dyu
      end if
!
      if (my_task == master_task) then
         read(25,rec=7,iostat=ioerr) TLAT_G
         close(25)
      endif

      call scatter_global(ANGLET, TLAT_G, master_task, distrb_clinic, &
                          field_loc_Center, field_type_angle)
      deallocate(TLAT_G,TLON_G)

      where (HTS <= c0) HTS = c1
      where (HTW <= c0) HTW = c1
      where (HUN <= c0) HUN = c1
      where (HUE <= c0) HUE = c1
      where (DXU <= c0) DXU = c1
      where (DYU <= c0) DYU = c1
      where (DXT <= c0) DXT = c1
      where (DYT <= c0) DYT = c1
      call POP_HaloUpdate(HTS, POP_haloClinic, POP_gridHorzLocSface, &
                          POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
      call POP_HaloUpdate(HTW, POP_haloClinic, POP_gridHorzLocWface, &
                          POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
      call POP_HaloUpdate(HUN, POP_haloClinic, POP_gridHorzLocWface, &
                          POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
      call POP_HaloUpdate(HUE, POP_haloClinic, POP_gridHorzLocSface, &
                          POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
      call POP_HaloUpdate(DXU, POP_haloClinic, POP_gridHorzLocSWcorner, &
                          POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
      call POP_HaloUpdate(DYU, POP_haloClinic, POP_gridHorzLocSWcorner, &
                          POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
      call POP_HaloUpdate(DXT, POP_haloClinic, POP_gridHorzLocCenter, &
                          POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
      call POP_HaloUpdate(DYT, POP_haloClinic, POP_gridHorzLocCenter, &
                          POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
   end if

!-----------------------------------------------------------------------
!EOC

 end subroutine read_horiz_grid

!***********************************************************************
!BOP
! !IROUTINE: vert_grid_internal
! !INTERFACE:

!***********************************************************************
!BOP
! !IROUTINE: compute_dz
! !INTERFACE:

 subroutine compute_dz(depth,zlength,dz_sfc,dz_deep)

! !DESCRIPTION:
!  Computes a thickness profile and total depth given the
!  parameters for the thickness function
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), intent(in) :: &
      zlength,          &! gaussian parameter for thickness func
      dz_sfc,           &! thickness of surface layer
      dz_deep            ! thickness of deep ocean layers

! !OUTPUT PARAMETERS:

   real (r8), intent(out) :: &
      depth               ! depth based on integrated thicknesses

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: k

!-----------------------------------------------------------------------

   depth = c0

!  do k=1,km
!     dz(k) = dz_deep - (dz_deep - dz_sfc)*exp(-(depth/zlength)**2)
!     depth = depth + dz(k)
!  end do

!-----------------------------------------------------------------------
!EOC

 end subroutine compute_dz

!***********************************************************************
!BOP
! !IROUTINE: read_vert_grid
! !INTERFACE:

 subroutine read_vert_grid(vert_grid_file)

! !DESCRIPTION:
!  Reads in layer thicknesses from grid input file
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (char_len), intent(in) :: &
      vert_grid_file

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: &
      k,                 &! vertical level index
      nu,                &! i/o unit number
      ioerr               ! i/o error flag

!-----------------------------------------------------------------------
!
!  read vertical layer thickness from file
!
!-----------------------------------------------------------------------
!--------------------------------------------------------------
!     VERTICAL LAYERED PARAMETERS
!--------------------------------------------------------------
!     ZKP    DEPTHS OF BOX BOTTOMS ON MODEL GRID "T" BOXES (M)
!            (0 -25 -50 -75 ... -5600)
!     ZKT    DEPTHS OF BOX BOTTOMS ON MODEL GRID "W" BOXES (M)
!     DZP    D(ZKP)
!     ODZP   1/DZP
!     ODZT   1/DZT

   if (my_task == master_task) then
      open (25, file=trim(vert_grid_file),form='unformatted')
      read(25) zkp
      close(25)
   endif


      call broadcast_array(zkp, master_task)
      lev1=zkp
      DO K = 1,KM
         DZP (K) = ZKP (K) - ZKP (K +1)
         ZKT (K) = (ZKP (K) + ZKP (K +1))*0.5D0
         lev (k)= (ZKP (K) + ZKP (K +1))*0.5D0
         ODZP (K)= 1.0D0/ DZP (K)
      END DO

      ODZT (1)= 2.0D0* ODZP (1)
      DO K = 2,KM
         ODZT (K)= 1.0D0/ (ZKT (K -1) - ZKT (K))
      END DO

!-----------------------------------------------------------------------
!EOC

 end subroutine read_vert_grid

!***********************************************************************

 subroutine read_topography(topography_file,kmt_global)

! !DESCRIPTION:
!  Reads in KMT field from file
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (char_len), intent(in) :: &
      topography_file     ! input file containing KMT field

   logical(log_kind), intent(in) :: &
      kmt_global       ! flag for generating only global KMT field

#include <netcdf.inc> !for reading basin index filed

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: iblock
   integer (i4) :: &
      nu                ,&! i/o unit number
      ioerr             ,&! i/o error flag
      reclength           ! record length

!-----------------------------------------------------------------------
!
!  read global KMT field (for use in setting up domain)
!
!-----------------------------------------------------------------------

   if (kmt_global) then

      allocate(KMT_G(imt_global,jmt_global))

      INQUIRE(iolength=reclength) KMT_G
      if (my_task == master_task) then
         open(25, file=topography_file,status='old',form='unformatted', &
                  access='direct', recl=reclength, iostat=ioerr)
      endif

      call broadcast_scalar(ioerr, master_task)
      if (ioerr /= 0) call exit_LICOM(sigAbort, &
                                 'Error opening topography_file')

      if (my_task == master_task) then
         read(25, rec=1, iostat=ioerr) KMT_G
         close(25)
!        do j=1672,1676
!        do i=1813,1820
!            if (kmt_g(i,j) > 0 )  kmt_g(i,j) = max(kmt_g(i,j),20)
!        end do
!        end do
!        do j=1672,1678
!        do i=1835,1840
!            if (kmt_g(i,j) > 0 )  kmt_g(i,j) = max(kmt_g(i,j),37)
!        end do
!        end do
!        do j=1805,1813
!        do i=2960,2970
!            if (kmt_g(i,j) > 0 )  kmt_g(i,j) = max(kmt_g(i,j),10)
!        end do
!        end do
!        do j=220,226
!        do i=860, 869
!            if (kmt_g(i,j) > 0 )  kmt_g(i,j) = max(kmt_g(i,j),29)
!        end do
!        end do
!        do j=792, 796
!        do i=3080, 3090
!             kmt_g(i,j) = 54
!        end do
!        end do
!        kmt_g(3242,1937)= 25
!        kmt_g(3243,1937)= 25
!        kmt_g(3243,1938)= 25
!        kmt_g(3243,1939)= 25
!        kmt_g(3244,1939)= 25
!        kmt_g(3244,1940)= 25
!        kmt_g(3245,1940)= 25
!        kmt_g(3246,1940)= 25
!        kmt_g(3246,1941)= 25
!        kmt_g(3247,1941)= 25
!        kmt_g(3248,1941)= 25
!        kmt_g(3248,1942)= 25
!        kmt_g(3248,1943)= 25
!        kmt_g(3249,1943)= 25
!        kmt_g(3249,1944)= 25
!        kmt_g(3250,1944)= 25
!        kmt_g(3251,1944)= 25
!        kmt_g(3251,1945)= 25
!        kmt_g(3252,1945)= 25
!        kmt_g(3253,1945)= 25
      endif

      call broadcast_scalar(ioerr, master_task)
      if (ioerr /= 0) call exit_LICOM(sigAbort, &
                                 'Error reading topography_file')

      call broadcast_array(KMT_G, master_task)

!-----------------------------------------------------------------------
!
!  otherwise read KMT field from file
!
!-----------------------------------------------------------------------

   else
      call scatter_global(KMT, KMT_G, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
     do iblock = 1,  nblocks_clinic
        write(*,*) "mytid=", mytid, iblock, blocks_clinic(iblock)
     end do

      call boundary
   end if
!
!-----------------------------------------------------------------------
!EOC

 end subroutine read_topography


!-----------------------------------------------------------------------
!ZWP
 subroutine read_basin(basin_file)

! !DESCRIPTION:
!  Reads in the Basin index field from file
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

#include <netcdf.inc> !for reading basin index filed
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
  character (char_len), intent(in) :: basin_file     
  integer (i4) :: &
      ioerr             ,&! i/o error flag
      reclength           ! record length

!-----------------------------------------------------------------------
!
!  read basin index filed from NC file !zwp 2013-10-17
!
!-----------------------------------------------------------------------

      allocate(BASIN_G(imt_global,jmt_global))
!

      INQUIRE(iolength=reclength) BASIN_G
      if (my_task == master_task) then
         open(25, file=basin_file,status='old',form='unformatted', &
                  access='direct', recl=reclength, iostat=ioerr)
      endif

      call broadcast_scalar(ioerr, master_task)
      if (ioerr /= 0) call exit_LICOM(sigAbort, &
                                 'Error opening basin file')

      if (my_task == master_task) then
         read(25, rec=1, iostat=ioerr) BASIN_G
         close(25)
      endif
!
      call scatter_global(basin, basin_g, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!
!
 end subroutine read_basin
!ZWP
!***********************************************************************
!BOP
! !IROUTINE: landmasks
! !INTERFACE:

 subroutine landmasks

! !DESCRIPTION:
!  Calculates additional masks for land points at each depth level.
!  These include real masks for applying multiplicative masks
!  instead of logical masks and also KMT arrays for neighbor points.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  masks for surface ocean T or U points
!
!-----------------------------------------------------------------------

   where (KMT >= 1)
      CALCT = .true.
      RCALCT = c1
   elsewhere
      CALCT = .false.
      RCALCT = c0
   endwhere

   where (KMU >= 1)
      CALCU = .true.
      RCALCU = c1
   elsewhere
      CALCU = .false.
      RCALCU = c0
   endwhere

!-----------------------------------------------------------------------
!
!  depth level fields (KMT,KMU) to north,south,east,west
!
!-----------------------------------------------------------------------

   KMTN = eoshift(KMT,dim=2,shift=-1)
   KMTS = eoshift(KMT,dim=2,shift=+1)
   KMTE = eoshift(KMT,dim=1,shift=+1)
   KMTW = eoshift(KMT,dim=1,shift=-1)

   KMUN = eoshift(KMU,dim=2,shift=-1)
   KMUS = eoshift(KMU,dim=2,shift=+1)
   KMUE = eoshift(KMU,dim=1,shift=+1)
   KMUW = eoshift(KMU,dim=1,shift=-1)

   KMTEE = eoshift(KMT,dim=1,shift=2)
   KMTNN = eoshift(KMT,dim=2,shift=-2)

!-----------------------------------------------------------------------
!EOC

 end subroutine landmasks

!***********************************************************************
!BOP
! !IROUTINE: area_masks
! !INTERFACE:

 subroutine area_masks(mask_filename,info_filename)

! !DESCRIPTION:
!   This subroutine reads in file with regional area mask and
!   marginal seas defined
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character(*), intent(in) :: &
      mask_filename           ,&! name of file containing region masks
      info_filename             ! name of file containing region names

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: &
      k, n,              &! loop counters
      nu,                &! i/o unit number
      reclength,         &! record length of file
      ioerr,             &! i/o error flag
      region              ! region counter

   real (r8) :: &
      tmp_vol,  &! temporary volume
      sea_area, &! total volume of a particular sea
      sea_vol    ! total volume of a particular sea

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      WORK                ! temporary space

   integer (i4), dimension(:,:), allocatable :: &
      REGION_G            ! global-sized region mask

!-----------------------------------------------------------------------
!
!  read in regional area masks, including marginal seas.  then
!  calculate related variables.
!
!-----------------------------------------------------------------------

!  allocate(REGION_MASK(nx_block,ny_block,max_blocks_clinic))

!  call get_unit(nu)
!  if (my_task == master_task) then
!     allocate(REGION_G(imt_global,jmt_global))
!     inquire (iolength=reclength) REGION_G
!     open(nu, file=mask_filename,status='old',form='unformatted', &
!              access='direct', recl=reclength, iostat=ioerr)
!  endif

!  call broadcast_scalar(ioerr, master_task)
!  if (ioerr /= 0) call exit_LICOM(sigAbort, &
!                                'Error opening region mask file')

!  if (my_task == master_task) then
!     read(nu, rec=1, iostat=ioerr) REGION_G
!     close(nu)
!  endif
!  call release_unit(nu)

!  call broadcast_scalar(ioerr, master_task)
!  if (ioerr /= 0) call exit_LICOM(sigAbort, &
!                                'Error reading region mask file')

!  call scatter_global(REGION_MASK, REGION_G, master_task,distrb_clinic, &
!                      field_loc_center, field_type_scalar)
!  if (my_task == master_task) deallocate(REGION_G)

!  num_regions = global_maxval(abs(REGION_MASK), &
!                              distrb_clinic, field_loc_center)
!
!---------------------------------------------------------------------
!
!       open and read file which contains region names and marginal-
!         seas balancing information
!
!---------------------------------------------------------------------


!  if(info_filename == 'unknown_region_info') then
!    call exit_LICOM (sigAbort,'ERROR: unknown region_info_filename')
!  endif
!
!  region_info(:)%name   = 'unknown_region_name'
!  region_info(:)%area   = c0
!  region_info(:)%volume = c0
!
!  call get_unit(nu)
 
!  if (my_task == master_task) then
 
!    open(nu,file=info_filename,form='formatted', &
!         status='unknown', iostat=ioerr)
!
!    do n = 1,num_regions
!      read (nu,*) region_info(n)%number     , &
!                  region_info(n)%name       , &
!                  region_info(n)%ms_bal%lat , &
!                  region_info(n)%ms_bal%lon , &
!                  region_info(n)%ms_bal%area
!    enddo
!
!    close(nu)
 
!  endif 
!
!  call release_unit(nu)
!
!  call broadcast_scalar(ioerr, master_task)
!  if (ioerr /= 0) call exit_LICOM(sigAbort, &
!                                'Error reading region name file')

!  do n = 1,num_regions
!    call broadcast_scalar(region_info(n)%name       ,master_task)
!    call broadcast_scalar(region_info(n)%number     ,master_task)
!    call broadcast_scalar(region_info(n)%ms_bal%lat ,master_task)
!    call broadcast_scalar(region_info(n)%ms_bal%lon ,master_task)
!    call broadcast_scalar(region_info(n)%ms_bal%area,master_task)
!  enddo
 
 
!---------------------------------------------------------------------
!
!       determine if region is a marginal sea
!
!---------------------------------------------------------------------
!  num_ms = 0
!  do n = 1,num_regions
!    if (region_info(n)%number < 0 )  then
!        num_ms = num_ms + 1
!        region_info(n)%marginal_sea = .true.
!    else
!        region_info(n)%marginal_sea = .false.
!    endif
!  enddo
!
!  if (num_ms > max_ms) then
!    if (my_task == master_task) then
!      write(stdout,*)'area_masks: maximum number of marginal seas exceeded'   
!    endif
!    call exit_LICOM(sigAbort, &
!                 'ERROR: must increase max_ms in module grid.F')
!  endif


!-----------------------------------------------------------------------
!
!  a negative value of REGION_MASK designates a marginal sea.  if
!  a region is a marginal sea, calculate the area and volume
!
!-----------------------------------------------------------------------

!  area_t_marg = c0
!  volume_t_marg = c0
!  volume_t_marg_k = c0

!  do region = 1, num_regions

!     WORK = merge(TAREA, c0, 1 <= KMT .and. REGION_MASK == -region)
!     sea_area = global_sum(WORK, distrb_clinic, field_loc_center)
!     area_t_marg = area_t_marg + sea_area

!     sea_vol = c0
!     do k = 1, km
!        WORK = merge(TAREA*dz(k), c0, &
!                     k <= KMT .and. REGION_MASK == -region)
!        tmp_vol = global_sum(WORK, distrb_clinic, field_loc_center)
!        sea_vol = sea_vol + tmp_vol
!        volume_t_marg_k(k) = volume_t_marg_k(k) + tmp_vol
!     enddo
!     volume_t_marg = volume_t_marg + sea_vol

!     if(sea_area /= c0) then
!       if (.not. region_info(region)%marginal_sea) then
!         call exit_LICOM (sigAbort,'ERROR: marginal-sea mismatch')
!       endif
!       region_info(region)%area   = sea_area
!       region_info(region)%volume = sea_vol
!     endif

!     if (my_task == master_task .and. sea_area /= c0) then
!        write(stdout,"('Region #',i2,' is a marginal sea')") region
!        write(stdout, &
!            "('  area (km^2) = ',e12.5, '  volume (km^3) = ',e12.5)") &
!            sea_area*1.0e-10_r8, sea_vol*1.0e-15_r8
!     endif

!  enddo


!---------------------------------------------------------------------
!
!       document regional information
!
!---------------------------------------------------------------------
 
!  if (my_task == master_task) then
!    write(stdout,blank_fmt)
!    write(stdout,1003)
!
!    do n = 1,num_regions
!     if (region_info(n)%marginal_sea) then
!       write(stdout,1004) region_info(n)%number                   , &
!                          trim(region_info(n)%name)               , &
!                          region_info(n)%area  *1.0e-10_r8        , &
!                          region_info(n)%volume*1.0e-15_r8      
!     else
!       write(stdout,1004) region_info(n)%number, trim(region_info(n)%name)
!     endif
!    enddo
!
!    write(stdout,blank_fmt)
!    write(stdout,1005)
 
!    do n = 1,num_regions
!     if (region_info(n)%marginal_sea) then
!       write(stdout,1006) region_info(n)%number        , &
!                          trim(region_info(n)%name)    , &
!                          region_info(n)%ms_bal%lat    , &
!                          region_info(n)%ms_bal%lon    , &
!                          region_info(n)%ms_bal%area
!     endif
!    enddo
!
!  endif 
!
!  if (my_task == master_task) then
!    write(stdout,blank_fmt)
!  endif

!  do n = 1,num_regions
!    if (region_info(n)%marginal_sea) then
!        region_info(n)%area  =region_info(n)%area  *1.0e-4_r8
!        region_info(n)%volume=region_info(n)%volume*1.0e-6_r8
!    endif
!  enddo


1003  format (  30x, '+', 23('-'),'+'                                     , &
              /,30x, '|',  2x,  'Marginal Seas Only   |'                  , &
              /,2x,'Region', 8x, 'Region ',   7x, '|',  2x,'Area'         , &
              8x,'Volume   |'                                             , &
              /,2x,'Number', 8x,  'Name',10x, '| (km^2)', 7x, '(km^3)'    , &
              3x, '|'                                                     , &
              /, 2x,6('-'), 1x,19('-'),2x,11('-'),2x                      , &
              11('-') )
1004  format (2x, i4, a22, 2(1pe13.5) )
1005  format (/,3x, ' Marginal Sea (E+P+M+R) Balancing Information'       , &
              /,2x,'Region', 8x, 'Region ',  8x, 'Bal.'                   , &
              4x,'Bal.',  4x, 'Search'                                    , &
              /,2x,'Number', 8x,  'Name',11x, 'Lat'  , 5x, 'Lon'          , &
              5x, 'Area'                                                  , &
              /,47x, '(cm^2)'                                             , &
              /, 2x,6('-'), 1x,19('-'),2x, 6('-'),2x                      , &
               6('-'), 2x, 11('-')  )                      
1006  format (1x, i4, a20, 3x,2(0pf8.2),   1pe13.5  )
!-----------------------------------------------------------------------
!EOC

 end subroutine area_masks

!***********************************************************************

 subroutine cf_area_avg

!-----------------------------------------------------------------------
!
!  calculate the coefficients of the 4-point stencils for
!  area-weighted averaging at T and U points
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  calculate {central,s,w,sw} coefficients for area-averaging
!  to T points.
!
!-----------------------------------------------------------------------


   AU0  = p25
   AUS  = p25
   AUW  = p25
   AUSW = p25

!-----------------------------------------------------------------------
!
!  calculate {central,n,e,ne} coefficients for area-averaging
!   to U points.
!
!-----------------------------------------------------------------------

   AT0  = UAREA
   ATN  = eoshift(UAREA,dim=2,shift=-1)
   ATE  = eoshift(UAREA,dim=1,shift=+1)
   ATNE = eoshift(ATE  ,dim=2,shift=-1)

   AT0  = AT0 *p25*TAREA_R
   ATN  = ATN *p25*TAREA_R
   ATE  = ATE *p25*TAREA_R
   ATNE = ATNE*p25*TAREA_R
!YU
   if (trim(horiz_grid_opt) == 'lat_lon') then
      AT0 = P25
      ATN = P25
      ATE = P25
      ATNE= P25
   end if
!YU

!-----------------------------------------------------------------------

 end subroutine cf_area_avg

!***********************************************************************
!BOP
! !IROUTINE: calc_upoints
! !INTERFACE:

 subroutine calc_upoints(errorCode)

! !DESCRIPTION:
!  Calculates lat/lon coordinates of T points from U points
!  using a simple average of four neighbors in Cartesian 3d space.
!
! !REVISION HISTORY:
!  same as module
!
! !OUTPUT PARAMETERS:

   integer (i4), intent(out) :: &
      errorCode                ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: i,j,n

   real (r8) ::                   &
      xc,yc,zc,xs,ys,zs,xw,yw,zw, &! Cartesian coordinates for
      xsw,ysw,zsw,tx,ty,tz,da      !    nbr points

   type (block) :: &
      this_block    ! block info for this block

!-----------------------------------------------------------------------
!
!  TLAT, TLON are southwest 4-point averages of ULAT,ULON
!  for general grids, must drop into 3-d Cartesian space to prevent
!  problems near the pole
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(n,i,j,this_block,xsw,ysw,zsw,xw,yw,zw,xs,ys,zs,xc,yc,zc, &
   !$OMP                     tx,ty,tz,da)

   do n=1,nblocks_clinic
      this_block = get_block(blocks_clinic(n),n)

      do j=1,ny_block-1
      do i=2,nx_block

!        !***
!        !*** set up averaging weights
!        !***

         !wt0  = AT0 (i,j,n)
         !wts  = ATS (i,j,n)
         !wtw  = ATW (i,j,n)
         !wtsw = ATSW(i,j,n)

!        !***
!        !*** convert neighbor U-cell coordinates to 3-d Cartesian coordinates 
!        !*** to prevent problems with averaging near the pole
!        !***

         zsw = cos(TLAT(i-1,j+1,n))
         xsw = cos(TLON(i-1,j+1,n))*zsw
         ysw = sin(TLON(i-1,j+1,n))*zsw
         zsw = sin(TLAT(i-1,j+1,n))

         zs  = cos(TLAT(i  ,j+1,n))
         xs  = cos(TLON(i  ,j+1,n))*zs
         ys  = sin(TLON(i  ,j+1,n))*zs
         zs  = sin(TLAT(i  ,j+1,n))

         zw  = cos(TLAT(i-1,j  ,n))
         xw  = cos(TLON(i-1,j  ,n))*zw
         yw  = sin(TLON(i-1,j  ,n))*zw
         zw  = sin(TLAT(i-1,j  ,n))

         zc  = cos(TLAT(i  ,j  ,n))
         xc  = cos(TLON(i  ,j  ,n))*zc
         yc  = sin(TLON(i  ,j  ,n))*zc
         zc  = sin(TLAT(i  ,j  ,n))

!        !***
!        !*** straight 4-point average to T-cell Cartesian coords
!        !***

         tx = p25*(xc + xs + xw + xsw)
         ty = p25*(yc + ys + yw + ysw)
         tz = p25*(zc + zs + zw + zsw)

!        !***
!        !*** convert to lat/lon in radians
!        !***

         da = sqrt(tx**2 + ty**2 + tz**2)

         ULAT(i,j,n) = asin(tz/da)

         if (tx /= c0 .or. ty /= c0) then
            ULON(i,j,n) = atan2(ty,tx)
         else
            ULON(i,j,n) = c0
         endif
           
      end do
      end do

      !***
      !*** for bottom row of domain where sw 4pt average is not valid,
      !*** extrapolate from interior
      !*** NOTE: THIS ASSUMES A CLOSED NORTH BOUNDARY - WILL NOT
      !***       WORK CORRECTLY FOR CYCLIC OPTION
      !***

      if (this_block%j_glob(this_block%je) == jmt_global) then
         do i=this_block%ib,this_block%ie
            ULON(i,this_block%je,n) = ULON(i,this_block%je-1,n)
            ULAT(i,this_block%je,n) = c2*ULAT(i,this_block%je-1,n) - &
                                         ULAT(i,this_block%je-2,n)
         end do
      endif

      where (ULON(:,:,n) > pi2) ULON(:,:,n) = ULON(:,:,n) - pi2
      where (ULON(:,:,n) < c0 ) ULON(:,:,n) = ULON(:,:,n) + pi2

   end do
   !$OMP END PARALLEL DO
!
!
   if (trim(horiz_grid_opt) == 'lat_lon') then
      do n = 1, nblocks_clinic
      do j = 1, jmt-1
      do i = 2, imt-1
         ulat(i,j,n) = p5*(tlat(i,j,n)+tlat(i,j+1,n))
         ulon(i,j,n) = p5*(tlon(i-1,j,n)+tlon(i,j,n))
      end do
      end do
      end do
   end if
!
!-----------------------------------------------------------------------
!
!  Update boundaries
!
!-----------------------------------------------------------------------

   
   call POP_HaloUpdate(ULAT, POP_haloClinic, POP_gridHorzLocSWcorner, & 
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

   if (errorCode /= 0) then
      call LICOM_ErrorSet(errorCode, &
         'calc_upoints: error updating tlat halo')
      return
   endif

   call POP_HaloUpdate(ULON, POP_haloClinic, POP_gridHorzLocSWcorner, & 
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

   if (errorCode /= 0) then
      call LICOM_ErrorSet(errorCode, &
         'calc_upoints: error updating tlon halo')
      return
   endif

!-----------------------------------------------------------------------

 end subroutine calc_upoints

!***********************************************************************
!BOP
! !IROUTINE: ugrid_to_tgrid
! !INTERFACE:

 subroutine ugrid_to_tgrid(ARRAY_TGRID, ARRAY_UGRID, iblock,k)

! !DESCRIPTION:
!  Interpolates values at U points on a B grid to T points.
!  Note that ghost cells are not updated.
!  Also note that the input array is assumed to be in the baroclinic
!  distribution (where the stencil weights are defined).
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (i4), intent(in) :: &
     iblock,k                  ! index for block in baroclinic distrb

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
     ARRAY_UGRID             ! field on U points

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: &
     ARRAY_TGRID             ! field on T points

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: &
     i,j                   ! dummy indices
   real (r8) :: epsln

!-----------------------------------------------------------------------
!
!  southwest 4pt average
!
!-----------------------------------------------------------------------
   epsln =1.0D-25

   do j=2,ny_block
   do i=1,nx_block-1

     ARRAY_TGRID(i,j) = vit(i,j,k,iblock)*(AT0 (i,j,iblock)*ARRAY_UGRID(i  ,j  ) + &
                        ATN (i,j,iblock)*ARRAY_UGRID(i  ,j-1) + &
                        ATE (i,j,iblock)*ARRAY_UGRID(i+1,j  ) + &
                        ATNE(i,j,iblock)*ARRAY_UGRID(i+1,j-1))/  &
                        (viv(i,j,k,iblock)*at0(i,j,iblock)   +viv(i,j-1,k,iblock)*atn(i,j,iblock) + &
                         viv(i+1,j,k,iblock)*ate(i,j,iblock) +viv(i+1,j-1,k,iblock)*atne(i,j,iblock) +epsln)

   end do
   end do

   ARRAY_TGRID(:,1) = c0
   ARRAY_TGRID(nx_block,:) = c0

!-----------------------------------------------------------------------
!EOC

 end subroutine ugrid_to_tgrid

!***********************************************************************
!BOP
! !IROUTINE: tgrid_to_ugrid
! !INTERFACE:

 subroutine tgrid_to_ugrid(ARRAY_UGRID, ARRAY_TGRID, iblock)

! !DESCRIPTION:
!  Interpolates values at T points on a B grid to U points.
!  Note that ghost cells are not updated.
!  Also note that the input array is assumed to be in the baroclinic
!  distribution (where the stencil weights are defined).
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (i4), intent(in) :: &
     iblock                  ! index for block in baroclinic distrb

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
     ARRAY_TGRID    ! field on T points

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: &
     ARRAY_UGRID    ! field on U points

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: i,j                 ! dummy indices

!-----------------------------------------------------------------------
!
!  northeast 4pt average
!
!-----------------------------------------------------------------------


   do j=1,ny_block-1
   do i=2,nx_block

      ARRAY_UGRID(i,j) = AU0 (i,j,iblock)*ARRAY_TGRID(i  ,j  ) + &
                         AUS (i,j,iblock)*ARRAY_TGRID(i  ,j+1) + &
                         AUW (i,j,iblock)*ARRAY_TGRID(i-1,j  ) + &
                         AUSW(i,j,iblock)*ARRAY_TGRID(i-1,j+1)

    end do
    end do
    ARRAY_UGRID(:,ny_block) = c0
    ARRAY_UGRID(1,:) = c0

!-----------------------------------------------------------------------
!EOC

 end subroutine tgrid_to_ugrid

!***********************************************************************

 subroutine calc_coeff
!
      real(r8) eps, abcd
      integer :: iblock, ErrorCode
      type (block) :: this_block
!
!$OMP PARALLEL DO PRIVATE(iblock,j,i,abcd)
      do iblock = 1, nblocks_clinic
!
      do j = 1,jmt
      do i=  1,imt
         snlat(i,j,iblock)= SIGN (1.0D0, ulat(i,j,iblock))
      end do
      end do
!
      do j = 1,jmt
      do i=  1,imt
         EPS = 0.5D0* FCOR(I,J,IBLOCK)* DTC
         EPEA(I,J,IBLOCK) = 1.0D0/ (1.0D0+ EPS * EPS)
         EPEB(I,J,IBLOCK) = EPS / (1.0D0+ EPS * EPS)
         EPS = FCOR(I,J,IBLOCK)* DTC
         EPLA(I,J,IBLOCK) = 1.0D0/ (1.0D0+ EPS * EPS)
         EPLB(I,J,IBLOCK) = EPS / (1.0D0+ EPS * EPS)
         EPS = 0.5D0* FCOR(I,J,IBLOCK)* DTB
         EBEA(I,J,IBLOCK) = 1.0D0/ (1.0D0+ EPS * EPS)
         EBEB(I,J,IBLOCK) = EPS / (1.0D0+ EPS * EPS)
         EPS = FCOR(I,J,IBLOCK)* DTB
         EBLA(I,J,IBLOCK) = 1.0D0/ (1.0D0+ EPS * EPS)
         EBLB(I,J,IBLOCK) = EPS / (1.0D0+ EPS * EPS)
      end do
      end do
!
      DO J = 1,JMT
         DO I = 1,IMT
!
            ABCD = 0.0
            DO K = 1,KM
               ABCD = ABCD+ VIT (I,J,K,iblock)* DZP (K)
            END DO
            IF (ABCD > 0.0)THEN
               OHBT (I,J,iblock)= 1.0D0/ ABCD
            ELSE
               OHBT (I,J,iblock)= 0.0D0
            END IF
!
            ABCD = 0.0
            DO K = 1,KM
               ABCD = ABCD+ VIV (I,J,K,iblock)* DZP (K)
            END DO
            DZPH (I,J,iblock)= ABCD
            IF (ABCD > 0.0)THEN
               OHBU (I,J,iblock)= 1.0D0/ ABCD
            ELSE
               OHBU (I,J,iblock)= 0.0D0
            END IF
         END DO
      END DO
   end do

   hbx =c0
   hby =c0
!$OMP PARALLEL DO PRIVATE(iblock,j,i)
   do iblock =1, nblocks_clinic
   do j=2,ny_block-1
   do i=2,nx_block-1
       hbx(i,j,iblock) = DXUR(i,j,iblock)*p5*(ht(i  ,j+1,iblock) - ht(i-1,j,iblock) - &
                                              ht(i-1,j+1,iblock) + ht(i  ,j,iblock))
       hby(i,j,iblock) = DYUR(i,j,iblock)*p5*(ht(i  ,j+1,iblock) - ht(i-1,j,iblock) + &
                                              ht(i-1,j+1,iblock) - ht(i  ,j,iblock))
   end do
   end do
   end do
!
  call POP_HaloUpdate(hbx, POP_haloClinic, POP_gridHorzLocSwcorner, &
                               POP_fieldKindVector, errorCode,         &
                               fillValue = 0.0_r8)
  call POP_HaloUpdate(hby, POP_haloClinic, POP_gridHorzLocSwcorner, &
                               POP_fieldKindVector, errorCode,         &
                               fillValue = 0.0_r8)
!
!$OMP PARALLEL DO PRIVATE(iblock,j,i)
    do iblock =1, nblocks_clinic
       do j= 1, jmt
       do i= 1, imt
          if (abs(FCORT(i,j,iblock)) > 0.0) then
             rrd1(i,j,iblock) = 2.0_r8/abs(FCORT(i,j,iblock))
          else
             rrd1(i,j,iblock) = 100000.0_r8
          end if
!
          if (abs(FCOR(i,j,iblock)) > 0.0) then
             rrd2(i,j,iblock) = 2.0_r8/abs(FCOR(i,j,iblock))
          else
             rrd2(i,j,iblock) = 100000.0_r8
          end if
!
          rrd1(i,j,iblock) = max(rrd1(i,j,iblock), 15000.0_r8)
          rrd2(i,j,iblock) = max(rrd2(i,j,iblock), 15000.0_r8)
          rrd1(i,j,iblock) = min(rrd1(i,j,iblock),100000.0_r8)
          rrd2(i,j,iblock) = min(rrd2(i,j,iblock),100000.0_r8)
       end do
       end do
    end do

 end subroutine calc_coeff

 end module grid

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

