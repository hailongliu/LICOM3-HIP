!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module hmix_del4

!BOP
! !MODULE: hmix_del4

! !DESCRIPTION:
!  This module contains routines for computing biharmonic horizontal
!  diffusion of momentum and tracers.
!
! !REVISION HISTORY:
!  SVN:$Id: hmix_del4.F90 28439 2011-05-18 21:40:58Z njn01 $

! !USES:

   use precision_mod
   use param_mod
   use LICOM_Error_Mod

   use blocks
   use msg_mod
   use distribution
   use constant_mod, only : c0, c1, c2, p5, radius, pi
   use grid
   use domain
   use broadcast
   use operators
   use tracer_mod, only : atb

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_del4u,  &
             init_del4t,  &
             hdiffu_del4, &
             hdifft_del4

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  private module variables
!
!  operator coefficients:
!
!  DT{N,S,E,W} = {N,S,E,W} coefficients of 5-point stencil for the
!                Del**2 operator before b.c.s have been applied
!  DU{C,N,S,E,W} = central and {N,S,E,W} coefficients of 5-point
!                  stencil for the Del**2 operator acting at U points
!                  (without metric terms that mix U,V)
!  DM{C,N,S,E,W} = central and {N,S,E,W} coefficients of 5-point
!                  stencil for the metric terms that mix U,V
!  DUM = central coefficient for metric terms that do not mix U,V
!
!-----------------------------------------------------------------------

   real (r8), dimension (:,:,:), allocatable :: &
      DTN,DTS,DTE,DTW,                          &
      DUC,DUN,DUS,DUE,DUW,                      &
      DMC,DMN,DMS,DME,DMW,                      &
      DUM,                                      &
      AHF,              &! variable mixing factor for tracer   mixing
      AMF,              &! variable mixing factor for momentum mixing
      ratio_dxy         !changing mixing !YYQ add

   real (r8) ::         &
      ah,               &! horizontal tracer   mixing coefficient
      am                 ! horizontal momentum mixing coefficient

   logical (log_kind) :: &
      lvariable_hmixu,   &! spatially varying mixing coeffs
      lvariable_hmixt     ! spatially varying mixing coeffs

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_del4u
! !INTERFACE:

 subroutine init_del4u

! !DESCRIPTION:
!  This routine calculates the coefficients of the 5-point stencils for
!  the biharmonic operator acting on momentum fields and also
!  calculates coefficients for all diffusive metric terms. See the
!  description under hdiffu for the form of the operator.
!
! !REVISION HISTORY:
!  same as module
!
! !OUTPUT PARAMETER:

   integer (i4) :: errorCode             ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::   &
      i,j,                 &! dummy loop indices
      iblock,              &! block index
      nml_error             ! error flag for namelist

   real (r8), dimension (:,:), allocatable :: &
      KXU,KYU,             &! metric factors
      DXKX,DYKY,DXKY,DYKX, &! d{x,y}k{x,y}
      WORK1,WORK2           ! temporary work space

   logical (log_kind) ::   &
      lauto_hmix,          &! flag to internally compute mixing coeff
      lvariable_hmix        ! flag to enable spatially varying mixing


   real (r8) ::            &
      amfmin, amfmax         ! min max mixing for varible mixing

!-----------------------------------------------------------------------
!
!  read input namelist to set options
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
!  set spatially variable mixing arrays if requested
!
!  this applies only to momentum or tracer mixing
!  with del2 or del4 options
!
!  functions describing spatial dependence of horizontal mixing
!  coefficients:   am -> am*AHM
!
!  for standard laplacian  mixing they scale like (cell area)**0.5
!  (note: these forms assume a global mesh with a grid line along
!  the equator, and the mixing coefficients are the set relative
!  to the average grid spacing at the equator - for cells of
!  this size AMF = AHF = 1.0).
!
!-----------------------------------------------------------------------
      allocate(AMF(nx_block,ny_block,nblocks_clinic))
!add 20200221
      allocate(ratio_dxy(nx_block,ny_block,nblocks_clinic))
!add 20200221
      do iblock = 1, nblocks_clinic
         AMF(:,:,iblock) = (uarea(:,:,iblock)/   &
                          ((2.0D0*radius*PI/DBLE(IMT_GLOBAL))**2.0))**1.5
!add 20200221
         ratio_dxy(:,:,iblock) = sqrt(uarea(:,:,iblock))/   &
                          (2.0D0*radius*PI/DBLE(IMT_GLOBAL))
!add 20200221
      end do
      AM  = AM_BIHAR
!-----------------------------------------------------------------------
!
!  calculate operator weights
!
!-----------------------------------------------------------------------

   allocate(DUC(nx_block,ny_block,nblocks_clinic), &
            DUN(nx_block,ny_block,nblocks_clinic), &
            DUS(nx_block,ny_block,nblocks_clinic), &
            DUE(nx_block,ny_block,nblocks_clinic), &
            DUW(nx_block,ny_block,nblocks_clinic), &
            DMC(nx_block,ny_block,nblocks_clinic), &
            DMN(nx_block,ny_block,nblocks_clinic), &
            DMS(nx_block,ny_block,nblocks_clinic), &
            DME(nx_block,ny_block,nblocks_clinic), &
            DMW(nx_block,ny_block,nblocks_clinic), &
            DUM(nx_block,ny_block,nblocks_clinic))

   allocate(KXU   (nx_block,ny_block), &
            KYU   (nx_block,ny_block), &
            DXKX  (nx_block,ny_block), &
            DYKY  (nx_block,ny_block), &
            DXKY  (nx_block,ny_block), &
            DYKX  (nx_block,ny_block), &
            WORK1 (nx_block,ny_block), &
            WORK2 (nx_block,ny_block))

   do iblock=1,nblocks_clinic

!-----------------------------------------------------------------------
!
!     calculate central and {N,S,E,W} coefficients for
!     Del**2 (without metric terms) acting on momentum.
!
!-----------------------------------------------------------------------

      WORK1 = HUN(:,:,iblock)/HTW(:,:,iblock)

      DUN(:,:,iblock) = WORK1*UAREA_R(:,:,iblock)
      DUS(:,:,iblock) = eoshift(WORK1,dim=2,shift=1)*UAREA_R(:,:,iblock)
      !*** DUN invalid now in northernmost ghost cell

      WORK1 = HUE(:,:,iblock)/HTS(:,:,iblock)

      DUE(:,:,iblock) = WORK1*UAREA_R(:,:,iblock)
      DUW(:,:,iblock) = eoshift(WORK1,dim=1,shift=-1)*UAREA_R(:,:,iblock)
      !*** DUE invalid now in easternmost ghost cell

!-----------------------------------------------------------------------
!
!     coefficients for metric terms in Del**2(U)
!     and for metric advection terms (KXU,KYU)
!
!-----------------------------------------------------------------------

      KXU = (HUE(:,:,iblock) - eoshift(HUE(:,:,iblock),dim=1,shift=-1))*&
            UAREA_R(:,:,iblock)
      !*** KXU invalid in easternmost ghost cell
      KYU = (eoshift(HUN(:,:,iblock),dim=2,shift=1) - HUN(:,:,iblock))*&
            UAREA_R(:,:,iblock)
      !*** KYU invalid in northernmost ghost cell

      WORK1 = (eoshift(HTW(:,:,iblock),dim=1,shift=1)-HTW(:,:,iblock))* &
              TAREA_R(:,:,iblock)  ! KXT
      !*** WORK1 invalid in westernmost ghost cell
      WORK2 = p5*(WORK1 + eoshift(WORK1,dim=2,shift=1))
      !*** WORK2 invalid in both easternmost and
      !***                       westernmost ghost cell
      DXKX = (WORK2 - eoshift(WORK2,dim=1,shift=-1))*DXUR(:,:,iblock)
      !*** DXKX invalid in east/west/northernmost ghost cells

      WORK2 = p5*(WORK1 + eoshift(WORK1,dim=1,shift=-1))
      !*** WORK2 invalid in west/northernmost ghost cells
      DYKX = (eoshift(WORK2,dim=2,shift=1) - WORK2)*DYUR(:,:,iblock)
      !*** DYKX invalid in east/west/northernmost ghost cells

      WORK1 = (HTS(:,:,iblock) -                         &
               eoshift(HTS(:,:,iblock),dim=2,shift=-1))* &
              TAREA_R(:,:,iblock)  ! KYT
      !*** WORK1 invalid in southernmost ghost cell
      WORK2 = p5*(WORK1 + eoshift(WORK1,dim=1,shift=-1))
      !*** WORK2 invalid in north/southernmost ghost cells
      DYKY = (eoshift(WORK2,dim=2,shift=1) - WORK2)*DYUR(:,:,iblock)

      WORK2 = p5*(WORK1 + eoshift(WORK1,dim=2,shift=1))
      DXKY = (WORK2 - eoshift(WORK2,dim=1,shift=-1))*DXUR(:,:,iblock)

      DUM(:,:,iblock) = -(DXKX + DYKY + c2*(KXU**2 + KYU**2))
      DMC(:,:,iblock) = DXKY - DYKX

!-----------------------------------------------------------------------
!
!     calculate central and {N,S,E,W} coefficients for
!     metric mixing terms which mix U,V.
!
!-----------------------------------------------------------------------

      DME(:,:,iblock) =  c2*KYU/(HTS(:,:,iblock) + &
                                 eoshift(HTS(:,:,iblock),dim=1,shift=-1))
      DMS(:,:,iblock) = -c2*KXU/(HTW(:,:,iblock) + &
                                 eoshift(HTW(:,:,iblock),dim=2,shift=1))

   end do

   DUC = -(DUN + DUS + DUE + DUW)               ! scalar laplacian
   DMW = -DME
   DMN = -DMS

!-----------------------------------------------------------------------
!
!  free up memory
!
!-----------------------------------------------------------------------

   deallocate(KXU, KYU,               &
              DXKX, DYKY, DXKY, DYKX, &
              WORK1, WORK2)


!-----------------------------------------------------------------------
!EOC

 end subroutine init_del4u

!***********************************************************************
!BOP
! !IROUTINE: init_del4t
! !INTERFACE:

 subroutine init_del4t

! !DESCRIPTION:
!  This routine reads parameters for biharmonic tracer mixing and
!  calculates the coefficients of the 5-point stencils for
!  the biharmonic operator acting on tracer fields.  See the hdifft
!  routine for a description of the operator.
!
! !REVISION HISTORY:
!  same as module
!
! !OUTPUT PARAMETERS:

   integer (i4) :: errorCode         ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,               &! dummy loop indices
      iblock,            &! block index
      nml_error           ! error flag for namelist

   real (r8), dimension (:,:), allocatable :: &
      WORK1                ! temporary work space

   logical (log_kind) :: &
      lauto_hmix,        &! true to automatically determine mixing coeff
      lvariable_hmix      ! true for spatially varying mixing


   real (r8) ::          &
      ahfmin, ahfmax      ! min max mixing for varible mixing

!-----------------------------------------------------------------------
!
!  read input namelist to set options
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  set spatially variable mixing arrays if requested
!
!  this applies only to momentum or tracer mixing
!  with del2 or del4 options
!
!  functions describing spatial dependence of horizontal mixing
!  coefficients:   ah -> ah*AHF
!
!  for standard laplacian  mixing they scale like (cell area)**0.5
!  (note: these forms assume a global mesh with a grid line along
!  the equator, and the mixing coefficients are the set relative
!  to the average grid spacing at the equator - for cells of
!  this size AMF = AHF = 1.0).
!
!-----------------------------------------------------------------------
      allocate(AHF(nx_block,ny_block,nblocks_clinic))
      do iblock = 1 ,nblocks_clinic
         AHF(:,:,iblock) = (tarea(:,:,iblock)/   &
                          ((2.0D0*radius*PI/DBLE(IMT_GLOBAL))**2.0))**1.5
      end do
   AH  = AH_BIHAR
!-----------------------------------------------------------------------
!
!  calculate {N,S,E,W} coefficients for Del**2 acting on tracer
!  fields (for tracers, the central coefficient is calculated as
!  minus the sum of these after boundary conditions are applied).
!
!-----------------------------------------------------------------------

   allocate(DTN(nx_block,ny_block,nblocks_clinic),  &
            DTS(nx_block,ny_block,nblocks_clinic),  &
            DTE(nx_block,ny_block,nblocks_clinic),  &
            DTW(nx_block,ny_block,nblocks_clinic))

   allocate(WORK1 (nx_block,ny_block))

   do iblock=1,nblocks_clinic

      WORK1 = HTS(:,:,iblock)/HUE(:,:,iblock)

      DTS(:,:,iblock) = WORK1*TAREA_R(:,:,iblock)
      DTN(:,:,iblock) = eoshift(WORK1,dim=2,shift=-1)* &
                        TAREA_R(:,:,iblock)

      WORK1 = HTW(:,:,iblock)/HUN(:,:,iblock)

      DTW(:,:,iblock) = WORK1*TAREA_R(:,:,iblock)
      DTE(:,:,iblock) = eoshift(WORK1,dim=1,shift=1)* &
                        TAREA_R(:,:,iblock)

   end do

!-----------------------------------------------------------------------
!
!  free up memory
!
!-----------------------------------------------------------------------

   deallocate(WORK1)

!-----------------------------------------------------------------------
!EOC

 end subroutine init_del4t

!***********************************************************************
!BOP
! !IROUTINE: hdiffu_del4
! !INTERFACE:

 subroutine hdiffu_del4(k,HDUK,HDVK,UMIXK,VMIXK,this_block)

! !DESCRIPTION:
!  This routine computes the horizontial diffusion of momentum
!  using a biharmonic ($\nabla^4$) operator, where the biharmonic
!  operator is implemented using a repeated application of the
!  Laplacian operator:
!  \begin{equation}
!  D_H(u) = \nabla^2(A_M\nabla^2(u))
!  \end{equation}
!  where
!  \begin{eqnarray}
!  \nabla^2(u) & = &  \Delta_x\delta_x[\Delta_y\delta_x(u)]/AREA
!                   + \Delta_y\delta_y[\Delta_x\delta_y(u)]/AREA \\
!              &   &  - u*[dxkx - dyky + 2(k_x^2 + k_y^2)]
!                + 2hy\delta_x(v) - 2k_x\delta_y(v) \nonumber \\
!  \nabla^2(v) & = & \Delta_x\delta_x[\Delta_y\delta_x(v)]/AREA
!                  + \Delta_y\delta_y[\Delta_x\delta_y(v)]/AREA \\
!              &   &  - v*[dxkx - dyky + 2(k_x^2 + k_y^2)]
!                - 2hy\delta_x(u) + 2k_x\delta_y(u) \nonumber
!  \end{eqnarray}
!
!  Boundary conditions are not explicitly imposed on
!  $\nabla^2(u,v)$ since $u = v = 0$ on the boundaries.
!  In this version of the model, the  boundary conditions
!  on $\nabla^2$ acting on $\nabla^2(u,v)$ are also that
!  $\nabla^2(u,v)$ vanishes at land points.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k  ! depth level index

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      UMIXK,             &! U velocity at level k and mix time level
      VMIXK               ! V velocity at level k and mix time level

   type (block), intent(in) :: &
      this_block             ! block info for this sub block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: &
      HDUK,              &! Hdiff(Ub) at level k
      HDVK                ! Hdiff(Vb) at level k

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,               &! loop indices
      bid                 ! local address of current block

   real (r8), dimension (nx_block,ny_block) :: &
      D2UK,D2VK,        &! intermediate Del**2 results
      CC,               &! central 5-point weight
      CN,CS,CE,CW,      &! additional weights for partial bottom cells
      HDIFFCFL           ! local hdiff cfl number for diagnostics
   real(r8), dimension(nx_block,ny_block) :: div_out,curl,gradx1,grady1,gradx2,grady2
   real(r8) :: dxdy,fact_1, am_factor(nx_block,ny_block,nblocks_clinic)

!-----------------------------------------------------------------------
!
!  biharmonic mixing
!
!-----------------------------------------------------------------------

   bid = this_block%local_id
!-----------------------------------------------------------------------
!
!  calculate Del**2(U,V) without metric terms that mix U,V
!  add metric terms that mix U,V, and zero fields at land points
!
!-----------------------------------------------------------------------
      am_factor = 1.0_r8
      call div(k,DIV_OUT,umixk,vmixk,this_block)
      call zcurl(k,CURL,umixk,vmixk,this_block)
      call grad(k, GRADX1, GRADY1, curl, this_block)
      call grad(k, GRADX2, GRADY2, div_out, this_block)
      do j=this_block%jb-1,this_block%je+1
      do i=this_block%ib-1,this_block%ie+1
          !dxdy = (sqrt(uarea(i,j,bid)))**5*20.0_r8
          dxdy = (sqrt(uarea(i,j,bid)))**5*45.0_r8
          am_factor(i,j,bid) = sqrt(gradx1(i,j)**2+gradx2(i,j)**2+  & 
                          grady1(i,j)**2+grady2(i,j)**2)*dxdy/abs(am*amf(i,j,bid))
      end do
      end do
!      if ( mod(mytid, 10) == 0 )then
!            close(150+mytid/10)
!            stop
!      end if
      am_factor = min(40.0_r8, am_factor)
      am_factor = max( 1.0_r8, am_factor)

   !*** add metric contribution to central coefficient
   CC = DUC(:,:,bid) + DUM(:,:,bid)
   D2UK = c0
   D2VK = c0

      do j=this_block%jb-1,this_block%je+1
      do i=this_block%ib-1,this_block%ie+1

         D2UK(i,j)  =(CC (i,j    )*UMIXK(i  ,j  ) +        &
                      DUN(i,j,bid)*UMIXK(i  ,j-1) +        &
                      DUS(i,j,bid)*UMIXK(i  ,j+1) +        &
                      DUE(i,j,bid)*UMIXK(i+1,j  ) +        &
                      DUW(i,j,bid)*UMIXK(i-1,j  ))+        &
                     (DMC(i,j,bid)*VMIXK(i  ,j  ) +        &
                      DMN(i,j,bid)*VMIXK(i  ,j-1) +        &
                      DMS(i,j,bid)*VMIXK(i  ,j+1) +        &
                      DME(i,j,bid)*VMIXK(i+1,j  ) +        &
                      DMW(i,j,bid)*VMIXK(i-1,j  ))
      end do
      end do

!
      do j=this_block%jb-1,this_block%je+1
      do i=this_block%ib-1,this_block%ie+1
         D2VK(i,j)  =(CC (i,j    )*VMIXK(i  ,j  ) +        &
                      DUN(i,j,bid)*VMIXK(i  ,j-1) +        &
                      DUS(i,j,bid)*VMIXK(i  ,j+1) +        &
                      DUE(i,j,bid)*VMIXK(i+1,j  ) +        &
                      DUW(i,j,bid)*VMIXK(i-1,j  ))-        &
                     (DMC(i,j,bid)*UMIXK(i  ,j  ) +        &
                      DMN(i,j,bid)*UMIXK(i  ,j-1) +        &
                      DMS(i,j,bid)*UMIXK(i  ,j+1) +        &
                      DME(i,j,bid)*UMIXK(i+1,j  ) +        &
                      DMW(i,j,bid)*UMIXK(i-1,j  ))
      end do
      end do

      where (k <= KMU(:,:,bid))
!20200221 YYQ
         !D2UK = AMF(:,:,bid)*D2UK
         !D2VK = AMF(:,:,bid)*D2VK
         D2UK = AM_FACTOR(:,:,bid)*AMF(:,:,bid)*D2UK
         D2VK = AM_FACTOR(:,:,bid)*AMF(:,:,bid)*D2VK
!20200221 YYQ
      elsewhere
         D2UK = c0
         D2VK = c0
      end where


!-----------------------------------------------------------------------
!
!  calculate Del**2(Del**2(U,V)) without metric terms that mix U,V
!  add metric terms, and zero fields at land points
!
!-----------------------------------------------------------------------

   HDUK = c0
   HDVK = c0

      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         HDUK(i,j)  = am*((CC (i,j    )*D2UK(i  ,j  ) +        &
                           DUN(i,j,bid)*D2UK(i  ,j-1) +        &
                           DUS(i,j,bid)*D2UK(i  ,j+1) +        &
                           DUE(i,j,bid)*D2UK(i+1,j  ) +        &
                           DUW(i,j,bid)*D2UK(i-1,j  ))+        &
                          (DMC(i,j,bid)*D2VK(i  ,j  ) +        &
                           DMN(i,j,bid)*D2VK(i  ,j-1) +        &
                           DMS(i,j,bid)*D2VK(i  ,j+1) +        &
                           DME(i,j,bid)*D2VK(i+1,j  ) +        &
                           DMW(i,j,bid)*D2VK(i-1,j  )))
      end do
      end do
!  if (mytid ==716 ) then
!     write(115,*) am,bid,cc(99,69)
!     write(115,*) dun(99,69,1),dus(99,69,1),due(99,69,1), duw(99,69,1)
!     write(115,*) umixk(99,69),umixk(99,68),umixk(99,70),umixk(98,69),umixk(100,69)
!     write(115,*) d2uk(99,69),d2uk(99,68),d2uk(99,70),d2uk(98,69),d2uk(100,69)
!     write(115,*) hduk(99,69),hduk(99,68),hduk(99,70),hduk(98,69),hduk(100,69)
!     write(115,*) am*cc(99,69)*d2uk(99,69),am*dun(99,69,1)*d2uk(99,68)
!     write(115,*) am*dus(99,69,1)*d2uk(99,70),am*due(99,69,1)*d2uk(100,68),am*duw(99,69,1)*d2uk(98,69)
!  end if 
      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         HDVK(i,j)  = am*((CC (i,j    )*D2VK(i  ,j  ) +        &
                           DUN(i,j,bid)*D2VK(i  ,j-1) +        &
                           DUS(i,j,bid)*D2VK(i  ,j+1) +        &
                           DUE(i,j,bid)*D2VK(i+1,j  ) +        &
                           DUW(i,j,bid)*D2VK(i-1,j  ))-        &
                          (DMC(i,j,bid)*D2UK(i  ,j  ) +        &
                           DMN(i,j,bid)*D2UK(i  ,j-1) +        &
                           DMS(i,j,bid)*D2UK(i  ,j+1) +        &
                           DME(i,j,bid)*D2UK(i+1,j  ) +        &
                           DMW(i,j,bid)*D2UK(i-1,j  )))
      end do
      end do


   where (k > KMU(:,:,bid))
      HDUK = c0
      HDVK = c0
   endwhere


 end subroutine hdiffu_del4

!***********************************************************************
!BOP
! !IROUTINE: hdifft_del4
! !INTERFACE:

!2020 YYQ
 subroutine hdifft_del4(k,D2TK,HDTK,TMIX,this_block)
 !subroutine hdifft_del4(k,HDTK,TMIX,this_block)
!2020 YYQ

! !DESCRIPTION:
!  This routine computes the horizontial diffusion of tracers
!  using a biharmonic ($\nabla^4$) operator, where the biharmonic
!  operator is implemented using a repeated application of the
!  Laplacian operator:
!  \begin{equation}
!  D_H(\phi) = \nabla^2(A_H\nabla^2(\phi))
!  \end{equation}
!  with
!  \begin{equation}
!  \nabla^2(\phi) = \Delta_x\delta_x[\Delta_y\delta_x(\phi)]/AREA
!                 + \Delta_y\delta_y[\Delta_x\delta_y(\phi)]/AREA
!  \end{equation}
!
!   The boundary conditions are the same on both
!   applications of the $\nabla^2$ operator:  the gradient of both
!   $T$ and $\nabla^2(T)$ vanishes across lateral boundaries.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k  


   type (block), intent(in) :: &
      this_block             ! block info for this sub block

! !OUTPUT PARAMETERS:

   real (r8), intent(in), dimension(nx_block,ny_block) :: &
      TMIX
   real (r8), intent(out), dimension(nx_block,ny_block) :: &
      HDTK                ! HDIFF(T) for tracer n at level k

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,n, iblock,              &! dummy loop indices
      bid                 ! local block address

   real (r8), dimension(nx_block,ny_block) :: &
      D2TK,              &! intermediate Del**2 result
      CC,CN,CS,CE,CW,    &! five point stencil coefficients
      WORK,              &! temp array for tavg quantities
      HDIFFCFL            ! local hdiff cfl number for diagnostics

!-----------------------------------------------------------------------
!
!  biharmonic mixing
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

!-----------------------------------------------------------------------
!
!  implement boundary conditions by setting
!  stencil coefficients to zero at land points.
!
!-----------------------------------------------------------------------

      CN = merge(DTN(:,:,bid), c0, (k <= KMTN(:,:,bid)) .and.  &
                                   (k <= KMT (:,:,bid)))
      CS = merge(DTS(:,:,bid), c0, (k <= KMTS(:,:,bid)) .and.  &
                                   (k <= KMT (:,:,bid)))
      CE = merge(DTE(:,:,bid), c0, (k <= KMTE(:,:,bid)) .and.  &
                                   (k <= KMT (:,:,bid)))
      CW = merge(DTW(:,:,bid), c0, (k <= KMTW(:,:,bid)) .and.  &
                                   (k <= KMT (:,:,bid)))

   CC = -(CN + CS + CE + CW)  ! central coefficient

   WORK = c0

!-----------------------------------------------------------------------
!
!     calculate Del**2(T)
!
!-----------------------------------------------------------------------


         do j=this_block%jb-1,this_block%je+1
         do i=this_block%ib-1,this_block%ie+1

            D2TK(i,j) = AHF(i,j,bid)*                    &
                       (CC(i,j)*TMIX(i  ,j  ) +      &
                        CN(i,j)*TMIX(i  ,j-1) +      &
                        CS(i,j)*TMIX(i  ,j+1) +      &
                        CE(i,j)*TMIX(i+1,j  ) +      &
                        CW(i,j)*TMIX(i-1,j  ))

         end do
         end do

!-----------------------------------------------------------------------
!
!     calculate Del**2(Del**2(T))
!     multiply by diffusivity
!
!-----------------------------------------------------------------------

      HDTK(:,:) = c0

      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie

         HDTK(i,j) = ah*(CC(i,j)*D2TK(i  ,j  ) +     &
                           CN(i,j)*D2TK(i  ,j-1) +     &
                           CS(i,j)*D2TK(i  ,j+1) +     &
                           CE(i,j)*D2TK(i+1,j  ) +     &
                           CW(i,j)*D2TK(i-1,j  ))

      end do
      end do
!    

!-----------------------------------------------------------------------
!EOC

 end subroutine hdifft_del4

!***********************************************************************

 end module hmix_del4

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
