!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module operators

!BOP
! !MODULE: operators
!
! !DESCRIPTION:
!  This module contains routines for various common
!  mathematical operators, including gradient, divergence, curl,
!  and a vertical velocity calculation.  Note that these routines
!  do {\em not} update ghost cells so results contain invalid
!  entries in some ghost cells.
!
! !REVISION HISTORY:
!  SVN:$Id: operators.F90 808 2006-04-28 17:06:38Z njn01 $

! !USES:

   use precision_mod
   use blocks
   use param_mod
   use domain
   use constant_mod
   use grid

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:
   public :: div,   &
             grad,  &
             zcurl

!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: div
! !INTERFACE:

 subroutine div(k,DIV_OUT,UX,UY,this_block)

! !DESCRIPTION:
!  This routine returns the divergence at T points (times the cell
!  area) of a vector field defined at U points.  The divergence
!  operator is defined as:
!
!  \begin{equation}
!  \nabla\cdot{\bf\rm u} = {1\over{\Delta_y}} \delta_x (\Delta_y u_x) +
!                          {1\over{\Delta_x}} \delta_y (\Delta_x u_y)
!  \end{equation}
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS

   integer (int_kind), intent(in) :: k   ! vertical level

   real (r8), dimension(nx_block,ny_block), intent(in) :: & 
      UX,UY              ! vector field defined at U-points

   type (block), intent(in) :: &
      this_block          ! block information for current block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: & 
      DIV_OUT            ! divergence at T-points times cell area

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,               &! dummy counters
      bid                 ! local block id

!-----------------------------------------------------------------------
!
!  compute divergence using a 4 point stencil
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

   DIV_OUT = c0

   do j=2,ny_block
   do i=1,nx_block-1
      if (k <= KMT(i,j,bid)) then
         DIV_OUT(i,j) = p5*((UX(i+1,j  )+ UX(i+1,j-1))*HTW(i+1,j,bid) -  &
                            (UX(i  ,j  )+ UX(i  ,j-1))*HTW(i  ,j,bid) +  &
                            (UY(i+1,j  )+ UY(i  ,j  ))*HTS(i  ,j,bid) -  &
                            (UY(i+1,j-1)+ UY(i  ,j-1))*HTS(i  ,j-1,bid))*TAREA_R(i,j,bid)
      endif
   end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine div

!***********************************************************************
!BOP
! !IROUTINE: grad
! !INTERFACE:

 subroutine grad(k, GRADX, GRADY, F, this_block)

! !DESCRIPTION:
!  This routine computes the gradient in i,j directions at U points
!  based on field defined at T-points.
!
!  \begin{eqnarray}
!     \nabla_x(F) &=& \delta_x F \\
!     \nabla_y(F) &=& \delta_y F
!  \end{eqnarray}
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k  ! vertical level

   real (r8), dimension(nx_block,ny_block), intent(in) :: & 
      F                  ! field defined at T points

   type (block), intent(in) :: &
      this_block          ! block information for current block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: & 
      GRADX,GRADY        ! gradient in (i,j) direction at U points

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,               &! dummy counters
      bid                 ! local block id

!-----------------------------------------------------------------------
!
!  compute gradient with a 4 point stencil
!
!-----------------------------------------------------------------------
      
   bid = this_block%local_id

   GRADX = c0
   GRADY = c0

   do j=1,ny_block-1
   do i=2,nx_block  
      if (k <= KMU(i,j,bid)) then
         GRADX(i,j) = DXUR(i,j,bid)*p5*(F(i  ,j+1) - F(i-1,j) - & 
                                        F(i-1,j+1) + F(i  ,j))
         GRADY(i,j) = DYUR(i,j,bid)*p5*(F(i  ,j+1) - F(i-1,j) + &
                                        F(i-1,j+1) - F(i  ,j))
      endif
   end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine grad

!***********************************************************************
!BOP
! !IROUTINE: zcurl
! !INTERFACE:

 subroutine zcurl(k,CURL,UX,UY,this_block)

! !DESCRIPTION:
!  This function returns the z-component of the curl of a vector
!  field defined at U points. 
!
!  \begin{equation}
!     \hat{\bf\rm z}\cdot\nabla\times{\bf\rm u} =
!     {1\over{\Delta_y}} \delta_x(\Delta_y u_y) -
!     {1\over{\Delta_x}} \delta_y(\Delta_x u_x)
!  \end{equation}
!
!  The result is actually multiplied by cell area and returned at 
!  T points. 
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k  ! vertical level

   real (r8), dimension(nx_block,ny_block), intent(in) :: & 
      UX,UY              ! vector field defined at U-points

   type (block), intent(in) :: &
      this_block          ! block information for current block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: & 
      CURL              ! z.curl(Ux,Uy) at T-points times cell area

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,               &! dummy counters
      bid                 ! local block index

!-----------------------------------------------------------------------
!
!  compute curl using stencil similar to divergence 4 point stencil
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

   CURL = c0

   do j=2,ny_block
   do i=2,nx_block
      if (k <= KMT(i,j,bid)) then
         CURL(i,j) = p5*(UY(i  ,j  )*DYU(i  ,j  ,bid) +  & 
                         UY(i  ,j-1)*DYU(i  ,j-1,bid) -  &
                         UY(i-1,j  )*DYU(i-1,j  ,bid) -  &
                         UY(i-1,j-1)*DYU(i-1,j-1,bid) -  &
                         UX(i  ,j  )*DXU(i  ,j  ,bid) -  &
                         UX(i-1,j  )*DXU(i-1,j  ,bid) +  &
                         UX(i  ,j-1)*DXU(i  ,j-1,bid) +  &
                         UX(i-1,j-1)*DXU(i-1,j-1,bid))
      endif
   end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine zcurl


 end module operators

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
