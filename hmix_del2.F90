!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module hmix_del2

!BOP
! !MODULE: hmix_del2

! !DESCRIPTION:
!  This module contains routines for computing Laplacian horizontal
!  diffusion of momentum and tracers.
!
! !REVISION HISTORY:
!  SVN:$Id: hmix_del2.F90 28439 2011-05-18 21:40:58Z njn01 $

! !USES:

   use precision_mod
   use LICOM_Error_Mod
   use pconst_mod, only : vit, viv
   use param_mod
   use blocks
   use msg_mod
   use distribution
   use domain
   use broadcast
   use constant_mod, only : c0, c1, c2, p5, radius, pi, degtorad
   use grid
   use tracer_mod, only : atb

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_del2u,  &
             init_del2t,  &
             hdiffu_del2, &
             hdifft_del2

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
      AMF                ! variable mixing factor for momentum mixing

   real (r8) ::        &
      ah,              &! horizontal tracer   mixing coefficient
      am                ! horizontal momentum mixing coefficient

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_del2u
! !INTERFACE:

 subroutine init_del2u

! !DESCRIPTION:
!  This routine calculates the coefficients of the 5-point stencils for
!  the $\nabla^2$ operator acting on momentum fields and also
!  calculates coefficients for all diffusive metric terms. See the
!  description under hdiffu for the form of the operator.
!
! !REVISION HISTORY:
!  same as module
!
! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::    &
      i,j,                  &! dummy loop indices
      iblock,               &! block index
      nml_error              ! error flag for namelist

   real (r8), dimension (:,:), allocatable :: &
      KXU,KYU,              &! metric factors
      DXKX,DYKY,DXKY,DYKX,  &! d{x,y}k{x,y}
      WORK1,WORK2            ! temporary work space

   real (r8) ::             &
      amfmin, amfmax         ! min max mixing for varible mixing


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
      do iblock = 1 ,nblocks_clinic
         AMF(:,:,iblock) = sqrt(uarea(:,:,iblock))/   &
                          (2.0D0*radius*pi/DBLE(IMT_GLOBAL))
!        do j =1, jmt
!        do i =1, imt
!           if ( ulat(i,j,iblock) < -15.0D0*degtorad) AMF(i,j,IBLOCK)= 1.0D0
!        end do
!        end do
      end do              

      AM  = AM_HOR

!-----------------------------------------------------------------------
!
!  calculate operator weights
!
!-----------------------------------------------------------------------

   allocate(DUC(nx_block,ny_block,nblocks_clinic),  &
            DUN(nx_block,ny_block,nblocks_clinic),  &
            DUS(nx_block,ny_block,nblocks_clinic),  &
            DUE(nx_block,ny_block,nblocks_clinic),  &
            DUW(nx_block,ny_block,nblocks_clinic),  &
            DMC(nx_block,ny_block,nblocks_clinic),  &
            DMN(nx_block,ny_block,nblocks_clinic),  &
            DMS(nx_block,ny_block,nblocks_clinic),  &
            DME(nx_block,ny_block,nblocks_clinic),  &
            DMW(nx_block,ny_block,nblocks_clinic),  &
            DUM(nx_block,ny_block,nblocks_clinic))

   allocate(KXU   (nx_block,ny_block),  &
            KYU   (nx_block,ny_block),  &
            DXKX  (nx_block,ny_block),  &
            DYKY  (nx_block,ny_block),  &
            DXKY  (nx_block,ny_block),  &
            DYKX  (nx_block,ny_block),  &
            WORK1 (nx_block,ny_block),  &
            WORK2 (nx_block,ny_block))

   do iblock=1,nblocks_clinic

!-----------------------------------------------------------------------
!
!     calculate central and {N,S,E,W} coefficients for
!     Del**2 (without metric terms) acting on momentum.
!
!-----------------------------------------------------------------------

      WORK1 = (HUN(:,:,iblock)/HTW(:,:,iblock))*p5*(AMF(:,:,iblock) + &
                            eoshift(AMF(:,:,iblock),dim=2,shift=-1))

      DUN(:,:,iblock) = WORK1*UAREA_R(:,:,iblock)
      DUS(:,:,iblock) = eoshift(WORK1,dim=2,shift=1)*UAREA_R(:,:,iblock)

      WORK1 = (HUE(:,:,iblock)/HTS(:,:,iblock))*p5*(AMF(:,:,iblock) + &
                            eoshift(AMF(:,:,iblock),dim=1,shift=1))

      DUE(:,:,iblock) = WORK1*UAREA_R(:,:,iblock)
      DUW(:,:,iblock) = eoshift(WORK1,dim=1,shift=-1)*UAREA_R(:,:,iblock)

!-----------------------------------------------------------------------
!
!     coefficients for metric terms in Del**2(U)
!     and for metric advection terms (KXU,KYU)
!
!-----------------------------------------------------------------------

      KXU = (HUE(:,:,iblock) - eoshift(HUE(:,:,iblock),dim=1,shift=-1))*&
            UAREA_R(:,:,iblock)
      KYU = (eoshift(HUN(:,:,iblock),dim=2,shift=1) - HUN(:,:,iblock))*&
            UAREA_R(:,:,iblock)

      WORK1 = (eoshift(HTW(:,:,iblock),dim=1,shift=1)-HTW(:,:,iblock))* &
              TAREA_R(:,:,iblock)  ! KXT

      WORK2 = p5*(WORK1 + eoshift(WORK1,dim=2,shift=1))*    &
              p5*(eoshift(AMF(:,:,iblock),dim=1,shift=1) + &
                  AMF(:,:,iblock))

      DXKX = (WORK2 - eoshift(WORK2,dim=1,shift=-1))*DXUR(:,:,iblock)

      WORK2 = p5*(WORK1 + eoshift(WORK1,dim=1,shift=-1))*    &
              p5*(eoshift(AMF(:,:,iblock),dim=2,shift=-1) + &
                  AMF(:,:,iblock))

      DYKX = (eoshift(WORK2,dim=2,shift=1) - WORK2)*DYUR(:,:,iblock)

      WORK1 = (HTS(:,:,iblock) -                         &
               eoshift(HTS(:,:,iblock),dim=2,shift=-1))* &
              TAREA_R(:,:,iblock)  ! KYT
      WORK2 = p5*(WORK1 + eoshift(WORK1,dim=1,shift=-1))*    &
              p5*(eoshift(AMF(:,:,iblock),dim=2,shift=-1) + &
                  AMF(:,:,iblock))

      DYKY = (eoshift(WORK2,dim=2,shift=1) - WORK2)*DYUR(:,:,iblock)

      WORK2 = p5*(WORK1 + eoshift(WORK1,dim=2,shift=1))*    &
              p5*(eoshift(AMF(:,:,iblock),dim=1,shift=-1) + &
                  AMF(:,:,iblock))

      DXKY = (WORK2 - eoshift(WORK2,dim=1,shift=-1))*DXUR(:,:,iblock)

      DUM(:,:,iblock) = -(DXKX + DYKY + &
                        c2*AMF(:,:,iblock)*(KXU**2 + KYU**2))
      DMC(:,:,iblock) = DXKY - DYKX

!-----------------------------------------------------------------------
!
!     calculate central and {N,S,E,W} coefficients for
!     metric mixing terms which mix U,V.
!
!-----------------------------------------------------------------------

      WORK1 = (eoshift(AMF(:,:,iblock),dim=2,shift= 1) -    &
               eoshift(AMF(:,:,iblock),dim=2,shift=-1))/    &
              (HTW(:,:,iblock) + eoshift(HTW(:,:,iblock),dim=2,shift=1))

      DME(:,:,iblock) =  (c2*AMF(:,:,iblock)*KYU + WORK1)/  &
                         (HTS(:,:,iblock) +                 &
                          eoshift(HTS(:,:,iblock),dim=1,shift=-1))

      WORK1 = (eoshift(AMF(:,:,iblock),dim=1,shift= 1) -    &
               eoshift(AMF(:,:,iblock),dim=1,shift=-1))/    &
              (HTS(:,:,iblock) + eoshift(HTS(:,:,iblock),dim=1,shift=-1))

      DMS(:,:,iblock) = -(c2*AMF(:,:,iblock)*KXU + WORK1)/  &
                         (HTW(:,:,iblock) +                 &
                          eoshift(HTW(:,:,iblock),dim=2,shift=1))

   end do

   DUC = -(DUN + DUS + DUE + DUW)               ! scalar laplacian
   DMW = -DME
   DMN = -DMS
!
     if ( trim(horiz_grid_opt) == 'lat_lon') then
       DMS= c0
       DMN= c0
       do iblock = 1, nblocks_clinic
       do j=1, jmt
       do i=1, imt
          DUM(i,j,iblock)= (c1-tan(ulat(i,j,iblock))**2)/radius/radius 
          DME(i,j,iblock)= sin(ulat(i,j,iblock))/(radius*dxu(i,j,iblock)*cos(ulat(i,j,iblock)))
       end do
       end do
       end do
       DMW = -DME
     end if
!
!-----------------------------------------------------------------------
!
!  free up memory
!
!-----------------------------------------------------------------------

   deallocate(KXU, KYU,                &
              DXKX, DYKY, DXKY, DYKX,  &
              WORK1, WORK2)

!-----------------------------------------------------------------------
!EOC

 end subroutine init_del2u

!***********************************************************************
!BOP
! !IROUTINE: init_del2t
! !INTERFACE:

 subroutine init_del2t

! !DESCRIPTION:
!  This routine reads parameters for Laplaciang tracer mixing and
!  calculates the coefficients of the 5-point stencils for
!  the $\nabla^2$ operator acting on tracer fields.  See the hdifft
!  routine for a description of the operator.
!
! !REVISION HISTORY:
!  same as module
!
! !OUTPUT PARAMETERS:


!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      i,j,                &! dummy loop indices
      iblock,             &! block index
      nml_error            ! error flag for namelist

   real (r8), dimension (:,:), allocatable :: &
      WORK1,WORK2          ! temporary work space

   logical (log_kind) ::  &
      lauto_hmix,         &! true to automatically determine mixing coeff
      lvariable_hmix       ! true for spatially varying mixing

   namelist /hmix_del2t_nml/ lauto_hmix, lvariable_hmix, ah

   real (r8) ::           &
      ahfmin, ahfmax       ! min max mixing for varible mixing
      


      allocate(AHF(nx_block,ny_block,nblocks_clinic))
      do iblock = 1 ,nblocks_clinic
         AHF(:,:,iblock) = sqrt(tarea(:,:,iblock))/   &
                          (2.0D0*RADIUS*PI/DBLE(IMT_GLOBAL))
      end do              
      AH  = AH_HOR !LPF20140423
      !AH  = AH_TRO

!-----------------------------------------------------------------------
!
!  calculate {N,S,E,W} coefficients for Del**2 acting on tracer
!  fields (for tracers, the central coefficient is calculated as
!  minus the sum of these after boundary conditions are applied).
!
!-----------------------------------------------------------------------

   allocate(DTN(nx_block,ny_block,nblocks_clinic), &
            DTS(nx_block,ny_block,nblocks_clinic), &
            DTE(nx_block,ny_block,nblocks_clinic), &
            DTW(nx_block,ny_block,nblocks_clinic))

   allocate(WORK1 (nx_block,ny_block), &
            WORK2 (nx_block,ny_block))

   do iblock=1,nblocks_clinic
      WORK1 = (HTS(:,:,iblock)/HUE(:,:,iblock))*p5*(AHF(:,:,iblock) + &
                            eoshift(AHF(:,:,iblock),dim=2,shift=1))

      DTS(:,:,iblock) = WORK1*TAREA_R(:,:,iblock)
      DTN(:,:,iblock) = eoshift(WORK1,dim=2,shift=-1)* &
                        TAREA_R(:,:,iblock)

      WORK1 = (HTW(:,:,iblock)/HUN(:,:,iblock))*p5*(AHF(:,:,iblock) + &
                            eoshift(AHF(:,:,iblock),dim=1,shift=-1))

      DTW(:,:,iblock) = WORK1*TAREA_R(:,:,iblock)
      DTE(:,:,iblock) = eoshift(WORK1,dim=1,shift=1)* &
                        TAREA_R(:,:,iblock)

   end do


!-----------------------------------------------------------------------
!
!  free up memory
!
!-----------------------------------------------------------------------

   deallocate(WORK1, WORK2)

!-----------------------------------------------------------------------
!EOC

 end subroutine init_del2t

!***********************************************************************
!BOP
! !IROUTINE: hdiffu_del2
! !INTERFACE:

 subroutine hdiffu_del2(k, HDUK, HDVK, UMIXK, VMIXK, this_block)

! !DESCRIPTION:
!  This routine computes the horizontial diffusion of momentum
!  using the Laplacian diffusion operator given by:
!  \begin{eqnarray}
!     \nabla\cdot A_M \nabla u & = &
!           {1\over{\Delta_y}}\delta_x
!           \left(\overline{A_M}^x \Delta_y \delta_x u \right)
!         + {1\over{\Delta_x}}\delta_y
!           \left(\overline{A_M}^y \Delta_x \delta_y u \right)
!           \nonumber \\
!      &  &- u\left(\delta_x k_x + \delta_y k_y +
!           2(k_x^2 + k_y^2)\right) \nonumber \\
!      &  &+ 2k_y \delta_x v - 2k_x \delta_y v \\
!     \nabla\cdot A_M \nabla v &=&
!           {1\over{\Delta_y}}\delta_x
!           \left(\overline{A_M}^x \Delta_y \delta_x v \right)
!         + {1\over{\Delta_x}}\delta_y
!           \left(\overline{A_M}^y \Delta_x \delta_y v \right)
!           \nonumber \\
!      &  &- v\left(\delta_x k_x + \delta_y k_y +
!           2(k_x^2 + k_y^2)\right) \nonumber \\
!      &  &+ 2k_y \delta_x u - 2k_x \delta_y u
!  \end{eqnarray}
!  where
!  \begin{equation}
!     k_x = {1\over{\Delta_y}}\delta_x\Delta_y
!  \end{equation}
!  and
!  \begin{equation}
!     k_y = {1\over{\Delta_x}}\delta_y\Delta_x
!  \end{equation}
!
!  Note that boundary conditions are not explicitly imposed on
!  since $u = v = 0$ on the boundaries.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k  ! depth level index

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      UMIXK,             &! U at level k and mixing time level
      VMIXK               ! V at level k and mixing time level

   type (block), intent(in) :: &
      this_block          ! block information for this subblock

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: &
      HDUK,              &! Hdiff(Ub) at level k
      HDVK                ! Hdiff(Vb) at level k

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      i,j,                &! loop indices
      bid                  ! local block address

   real (r8) ::          &
      cc,                &! center fivept weight
      cn, cs, ce, cw      ! other weights for partial bottom cells

   real (r8), dimension(nx_block,ny_block) :: &
      UTMP, VTMP,        &! modified velocities to use with topostress
      HDIFFCFL            ! for cfl number diagnostics

!-----------------------------------------------------------------------
!
!  laplacian mixing
!
!  calculate Del**2(U,V) without metric terms that mix U,V
!  add metric terms that mix U,V
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

   HDUK = c0
   HDVK = c0

!-----------------------------------------------------------------------
!
!  handle four cases individually to avoid unnecessary copies
!  these are all basic five point stencil operators - the topostress
!  option requires operating on a modified velocity while the
!  partial bottom cell case modifies the weights.
!
!-----------------------------------------------------------------------
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie

            !*** add metric contrib to central coeff
            cc = DUC(i,j,bid) + DUM(i,j,bid)

            HDUK(i,j) = am*((cc          *UMIXK(i  ,j  ) +  &
                             DUN(i,j,bid)*UMIXK(i  ,j-1) +  &
                             DUS(i,j,bid)*UMIXK(i  ,j+1) +  &
                             DUE(i,j,bid)*UMIXK(i+1,j  ) +  &
                             DUW(i,j,bid)*UMIXK(i-1,j  ))+  &
                            (DMC(i,j,bid)*VMIXK(i  ,j  ) +  &
                             DMN(i,j,bid)*VMIXK(i  ,j-1) +  &
                             DMS(i,j,bid)*VMIXK(i  ,j+1) +  &
                             DME(i,j,bid)*VMIXK(i+1,j  ) +  &
                             DMW(i,j,bid)*VMIXK(i-1,j  )))*viv(i,j,k,bid)

            HDVK(i,j) = am*((cc          *VMIXK(i  ,j  ) +  &
                             DUN(i,j,bid)*VMIXK(i  ,j-1) +  &
                             DUS(i,j,bid)*VMIXK(i  ,j+1) +  &
                             DUE(i,j,bid)*VMIXK(i+1,j  ) +  &
                             DUW(i,j,bid)*VMIXK(i-1,j  ))-  &
                            (DMC(i,j,bid)*UMIXK(i  ,j  ) +  &
                             DMN(i,j,bid)*UMIXK(i  ,j-1) +  &
                             DMS(i,j,bid)*UMIXK(i  ,j+1) +  &
                             DME(i,j,bid)*UMIXK(i+1,j  ) +  &
                             DMW(i,j,bid)*UMIXK(i-1,j  )))*viv(i,j,k,bid)
!     if (mytid ==8 .and. i > 2 .and. i < 4 .and. j >26 .and. j < 29 ) then
!               write(220,*) i,j,k, this_block%i_glob(i), this_block%j_glob(j)
!               write(220,*) duc(i,j,1),dum(i,j,1),dun(i,j,1),dus(i,j,1),due(i,j,1),duw(i,j,1)
!               write(220,*) dme(i,j,1),dmw(i,j,1)
!               write(220,*) dmn(i,j,1),dms(i,j,1),dmc(i,j,1)
!               write(220,*) umixk(i,j),vmixk(i,j), am
!               write(220,*) am*cc*UMIXK(i  ,j  ), am* DUN(i,j,bid)*UMIXK(i ,j-1) , am*DUS(i,j,bid)*UMIXK(i ,j+1)
!               write(220,*) am*DUE(i,j,bid)*UMIXK(i+1,j),am*DUW(i,j,bid)*UMIXK(i-1,j), am*DMC(i,j,bid)*VMIXK(i ,j ) 
!               write(220,*) am*DMN(i,j,bid)*VMIXK(i,j-1), am*DMS(i,j,bid)*VMIXK(i,j+1),am*DME(i,j,bid)*VMIXK(i+1,j )
!               write(220,*) am*DMW(i,j,bid)*VMIXK(i-1,j  )
!               write(220,*) am*cc*VMIXK(i  ,j  ), am* DUN(i,j,bid)*VMIXK(i ,j-1) , am*DUS(i,j,bid)*VMIXK(i ,j+1)
!               write(220,*) am*DUE(i,j,bid)*VMIXK(i+1,j),am*DUW(i,j,bid)*VMIXK(i-1,j), am*DMC(i,j,bid)*UMIXK(i ,j ) 
!               write(220,*) am*DMN(i,j,bid)*UMIXK(i,j-1), am*DMS(i,j,bid)*UMIXK(i,j+1),am*DME(i,j,bid)*UMIXK(i+1,j )
!               write(220,*) am*DMW(i,j,bid)*UMIXK(i-1,j  )
!               write(220,*) hduk(i,j), hdvk(i,j)
!           end if

         end do
         end do

!        if (mytid ==8) close(220)
!-----------------------------------------------------------------------
!
!  zero fields at land points
!
!-----------------------------------------------------------------------

   where (k > KMU(:,:,bid))
      HDUK = c0
      HDVK = c0
   endwhere

!-----------------------------------------------------------------------
!EOC

 end subroutine hdiffu_del2

!***********************************************************************
!BOP
! !IROUTINE: hdifft_del2
! !INTERFACE:

 subroutine hdifft_del2(k,HDTK,TMIX,this_block)

! !DESCRIPTION:
!  This routine computes the horizontial diffusion of tracers
!  using the Laplacian operator given by:
!  \begin{equation}
!     \nabla\cdot A_M \nabla \phi =
!           {1\over{\Delta_y}}\delta_x
!           \left(\overline{A_H}^x \Delta_y \delta_x \phi \right)
!         + {1\over{\Delta_x}}\delta_y
!           \left(\overline{A_H}^y \Delta_x \delta_y \phi \right)
!  \end{equation}
!  with the boundary conditions of zero gradients of tracers.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k


   type (block), intent(in) :: &
      this_block          ! block information for this subblock

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(in) ::  &
      TMIX
   real (r8), dimension(nx_block,ny_block), intent(out) ::  &
      HDTK                ! HDIFF(T) for tracer n at level k

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      i,j,n, iblock,      &! dummy tracer index
      bid                  ! local block address

   real (r8), dimension(nx_block,ny_block) :: CC,CN,CS,CE,CW     ! coeff of 5pt stencil for Del**2

!-----------------------------------------------------------------------
!
!  laplacian mixing
!
!  implement boundary conditions by setting
!  stencil coefficients to zero at land points.
!
!-----------------------------------------------------------------------

   bid = this_block%local_id


      CN = merge(DTN(:,:,bid), c0, (k <= KMTN(:,:,bid)) .and. &
                                   (k <= KMT (:,:,bid)))
      CS = merge(DTS(:,:,bid), c0, (k <= KMTS(:,:,bid)) .and. &
                                   (k <= KMT (:,:,bid)))
      CE = merge(DTE(:,:,bid), c0, (k <= KMTE(:,:,bid)) .and. &
                                   (k <= KMT (:,:,bid)))
      CW = merge(DTW(:,:,bid), c0, (k <= KMTW(:,:,bid)) .and. &
                                   (k <= KMT (:,:,bid)))

   CC = -(CN + CS + CE + CW)  ! central coefficient

!-----------------------------------------------------------------------
!
!  calculate Del**2(T) for each tracer n
!
!-----------------------------------------------------------------------

   HDTK = c0

   do j=this_block%jb,this_block%je
   do i=this_block%ib,this_block%ie
      HDTK(i,j) = ah*(CC(i,j)*TMIX(i  ,j  ) + &
                      CN(i,j)*TMIX(i  ,j-1) + &
                      CS(i,j)*TMIX(i  ,j+1) + &
                      CE(i,j)*TMIX(i+1,j  ) + &
                      CW(i,j)*TMIX(i-1,j  ))
   enddo
   enddo

!-----------------------------------------------------------------------
!EOC

 end subroutine hdifft_del2

!***********************************************************************

 end module hmix_del2

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
