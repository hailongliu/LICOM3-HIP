 module LICOM_Error_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
   use precision_mod
   use param_mod
   use msg_mod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: exit_LICOM, LICOM_ErrorSet

! !DEFINED PARAMETERS:

   integer (int_kind), parameter, public :: &
      sigExit  =  0,    &! signal for normal exit
      sigAbort = -1      ! signal for aborting (exit due to error)

!EOP
   integer , parameter, public :: LICOM_stdout = 6
!BOC
!EOC
!***********************************************************************

 contains


 subroutine LICOM_ErrorSet(errorCode, errorMsg)

! !DESCRIPTION:
!  This routine sets an error code to POP\_Fail and adds a message to
!  the error log for later printing.
!
! !REVISION HISTORY:
!  same as module


   integer (i4), intent(out) :: errorCode              ! Error code to set to fail

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      errorMsg               ! message to add to error log for printing

!-----------------------------------------------------------------------
!
!  Add error message to error log
!     
!-----------------------------------------------------------------------
       write(6,*) "LICOM error:", errorMsg
       close(6)

       errorCode = -1
!-----------------------------------------------------------------------
!EOC

 end subroutine LICOM_ErrorSet



 subroutine exit_LICOM(exit_mode, exit_message)


! !PUBLIC MEMBER FUNCTIONS:

!
! !DEFINED PARAMETERS:

   integer (int_kind), intent(in) :: &
     exit_mode    ! method for exiting (normal exit or abort)

   character (*), intent(in) :: &
     exit_message ! message to print before stopping

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      ierr,               &! error flag
      local_unit

!-----------------------------------------------------------------------
!
!  print message - must use unit 6 in place of stdout to
!  prevent circular dependence with io module
!
!-----------------------------------------------------------------------

     local_unit = 6

#ifndef CCSMCOUPLED
   if (my_task == master_task) then
#endif

      select case(exit_mode)
      case(sigExit)
         write (local_unit,'(a14)') 'LICOM exiting...'
      case(sigAbort)
         write (local_unit,'(a15)') 'LICOM aborting...'
      case default
         write (local_unit,'(a37)') 'LICOM exiting with unknown exit mode...'
      end select

      write (local_unit,*) exit_message
      call flush(local_unit)
#ifndef CCSMCOUPLED
   endif
#endif

!-----------------------------------------------------------------------
!
!  exit or abort the message-passing environment if required
!
!-----------------------------------------------------------------------

   stop

!-----------------------------------------------------------------------
!EOC

 end subroutine exit_LICOM

   end  module LICOM_Error_mod

