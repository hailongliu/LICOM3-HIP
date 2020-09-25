!  CVS: $Id: param_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module param_mod
#include <def-undef.h>
use precision_mod
!-------------------------------------------------------------------------------
!
! Author: Yongqiang YU  ( 12 Nov, 2002)
!
!-------------------------------------------------------------------------------
!YU  01/02/2013
      integer, parameter:: LICOM_blockSizeX = BLCKX ! size of block in first  horizontal dimension
      integer, parameter:: LICOM_blockSizeY = BLCKY ! size of block in second horizontal dimension
      integer, parameter:: max_blocks_clinic = MXBLCKS, max_blocks_tropic = MXBLCKS    !   in each distribution
      integer , parameter, public :: nghost = 2       ! number of ghost cells around each block
      integer , parameter, public :: nx_block = licom_BlockSizeX + 2*nghost,   &!  x,y dir including ghost
                                     ny_block = licom_BlockSizeY + 2*nghost     !  cells
!YU
      integer,parameter:: jmt_global=NJMT  ! Number of the End Grid for Tracer in Latitude.
      integer,parameter:: jmm_global=jmt_global-1
      integer,parameter:: imt_global=NIMT    ! Number of Grid Points in Longitude
      integer,parameter:: km=NKM     ! Number of Grid Points in Vertical Direction
!
      integer,parameter:: num_overlap=2 ! Number of overlapping grids for subdomain.

!Nummber of grids in the each subdomain
      integer,parameter:: jst=1     ! Number of the Strating Grid for Tracer in Latitude.
      integer,parameter:: jsm=jst+1 ! Number of the Strating Grid for Momentum in Latitude.
!ZHW should Add 1
      integer,parameter:: jet= ny_block
      integer,parameter:: jem=jet-1 ! Number of the End Grid for Momentum in Latitude.
      integer,parameter:: jmt= ny_block
      integer,parameter:: imt= nx_block
      integer :: j_loop             ! Loop index of J cycle for the each subdomain.
!
      integer,parameter:: imm_global=imt_global-1,imm=imt-1,jmm=jmt-1,kmp1=km+1,kmm1=km-1
      integer,parameter:: ntra=2    ! Number of Tracers
      integer:: i,j,k,m,n,ierr, mytid
      integer::jj_start,jj_end ! No Overlaped j ,North to South
      integer:: my_task, master_task

      real(r8) :: am_hor ! Horizontal laplacian viscosity 
      real(r8) :: ah_hor ! Horizontal laplacian diffusity
      real(r8) :: am_bihar ! Horizontal biharmonic viscosity 
      real(r8) :: ah_bihar ! Horizontal biharmonic diffusity
      integer :: num_cpl ! couple frequency 

!lhl090729
      integer,parameter:: s_imt=640,s_jmt=320
!lhl090729
!YU  01/02/2013

end module param_mod
