!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP

 module POP_CommMod

! !MODULE: POP_CommMod
! !DESCRIPTION:
!  This module contains necessary routines and variables to support
!  other parallel communication modules in POP.  In particular, this
!  module contains communicators, tags, task ids and other necessary
!  information and the routines to set them up.  In addition, several
!  utility routines for setting up the communication environment are
!  included.
!
! !REVISION HISTORY:
!  SVN:$Id: POP_CommMod.F90 72 2007-10-05 21:01:46Z pwjones $
!  2006-07-10: Phil Jones
!              add new communication module with new naming conventions
!              also adds new CCSM coupler interface based on merge
!                 between POP 2.0 code and CCSM POP from Nancy Norton
!
! !USES:

   use precision_mod
   use msg_mod, only : mpi_comm_ocn

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public  :: POP_CommCreateCommunicator

   

! !PUBLIC DATA MEMBERS:

   integer (i4), public :: &
      POP_communicator,         &! MPI communicator for ocn comms
      POP_myTask,               &! MPI task number for this task
      POP_masterTask             ! MPI task number for master task

   integer (i4), parameter, public :: &
      POP_mpitagHalo           = 1,       &! MPI tags for various
      POP_mpitagRedist         = 1000      ! communication patterns

!EOP
!BOC
!EOC
!***********************************************************************

 contains

 subroutine POP_CommCreateCommunicator(newCommunicator, numProcs)

! !DESCRIPTION:
!  This routine creates a separate communicator for a subset of
!  processors under default ocean communicator.
!
! !REVISION HISTORY:
!  same as module

! !INCLUDES:

! !INPUT PARAMETERS:

   integer (i4), intent(in) :: &
      numProcs          ! num of procs in new distribution

! !OUTPUT PARAMETERS:

   integer (i4), intent(out) :: &
      newCommunicator   ! new communicator for this distribution

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: &
     MPI_GROUP_OCN,         &! group of processors assigned to ocn
     MPI_GROUP_NEW           ! group of processors assigned to new dist

   integer (i4) :: &
     ierr                    ! error flag for MPI comms

   integer (i4), dimension(3) :: &
     range                   ! range of tasks assigned to new dist
                             !  (assumed 0,num_procs-1)

!-----------------------------------------------------------------------
!
!  determine group of processes assigned to distribution
!
!-----------------------------------------------------------------------
   POP_communicator = mpi_comm_ocn

   call MPI_COMM_GROUP (POP_Communicator, MPI_GROUP_OCN, ierr)

   range(1) = 0
   range(2) = numProcs-1
   range(3) = 1

!-----------------------------------------------------------------------
!
!  create subroup and communicator for new distribution
!  note: MPI_COMM_CREATE must be called by all procs in POP_Communicator
!
!-----------------------------------------------------------------------

   call MPI_GROUP_RANGE_INCL(MPI_GROUP_OCN, 1, range, &
                             MPI_GROUP_NEW, ierr)

   call MPI_COMM_CREATE (POP_Communicator, MPI_GROUP_NEW,  &
                         newCommunicator, ierr)

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_CommCreateCommunicator

 end module POP_CommMod
