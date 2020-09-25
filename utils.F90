#include "def-undef.h"
module utils
! Purpose:
! This module is devised for report errors, show in-format message,
!   handle errors for the netcdf-f90 interface, et al.
!
! Record of revisions:
!    Date             Programmer                Description of change  
!==============    ===================   =============================
! 2011/12/22         WenYu Huang                Original Code
! ...
!
!Any comments please send to huangwenyu@mail.tsinghua.edu.cn
  use precision_mod
  use param_mod
  contains

    !---------------------------------------------------------------------
    !                       Subroutine handle_err
    !---------------------------------------------------------------------
    !---------------------------------------------------------------------
    subroutine check_single_value_real_r8(var_name, var)
        character(len = *), intent(in) :: var_name
        real(r8), intent(in) :: var
        write(*, *), var_name, " : ", var
    end subroutine check_single_value_real_r8
    !---------------------------------------------------------------------
    !                  End of Subroutine check_single_value
    !---------------------------------------------------------------------
    !---------------------------------------------------------------------
    !               Subroutine check_single_value_real_r8_1d
    !---------------------------------------------------------------------
    subroutine check_single_value_real_r8_1d(var_name,var,num)
        character(len = *), intent(in) :: var_name
        integer, intent(in) :: num
        real(r8), intent(in) :: var(num) 
        write(*, *), var_name, " : ", var
    end subroutine check_single_value_real_r8_1d
    !---------------------------------------------------------------------
    !           End of Subroutine check_single_value_real_r8_1d
    !---------------------------------------------------------------------

end module utils
