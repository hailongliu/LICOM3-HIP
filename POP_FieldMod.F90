!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_FieldMod

!BOP
! !MODULE: POP_FieldMod
!

   implicit none
   private
   save

! !PUBLIC TYPES:

   ! Full field data type


   character (7), parameter, public :: POP_fieldKindUnknown  = 'unknown'
   character (6), parameter, public :: POP_fieldKindScalar   = 'scalar'
   character (6), parameter, public :: POP_fieldKindVector   = 'vector'
   character (5), parameter, public :: POP_fieldKindAngle    = 'angle'
   character (8), parameter, public :: POP_fieldKindNoUpdate = 'noUpdate'

   !*** identifiers for commonly-used POP data types

!BOC

 end module POP_FieldMod
