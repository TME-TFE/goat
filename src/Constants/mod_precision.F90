!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! General overarching module for precision definition. 

module mod_precision

    ! Initialize
    !============
    ! The usual
    implicit none
    save
    public 

    ! Precision types
    integer, parameter       :: R4      = selected_real_kind(6)
    integer, parameter       :: R8      = selected_real_kind(14)
    integer, parameter       :: R16     = selected_real_kind(33)

    integer, parameter       :: C32     = 32
    integer, parameter       :: C64     = 64
    integer, parameter       :: C128    = 128

    integer, parameter       :: I4      = selected_int_kind(4)
    integer, parameter       :: I8      = selected_int_kind(8)
    integer, parameter       :: I16     = selected_int_kind(16)

    ! Interfaces
    interface IsRealEqualToInt
        module procedure IsRealEqualToIntR8I8
    end interface

contains 
    
    ! Comparison subroutines
    function IsRealEqualToIntR8I8(realval, intval, tolfacopt) result(isequal)

        ! Description
        !============
        ! This routine checks whether a real value is equal to an integer
        ! value up to a certain precision relative to the precision of
        ! the real value. Useful to check e.g. if 1.000001 corresponds to 1
        ! and deal with machine precision effects. 

        ! Declare variables
        !==================
        ! Arguments
        real(R8), intent(in)            :: realval 
        integer(I8), intent(in)         :: intval 
        real(R8), intent(in), optional  :: tolfacopt
        logical                         :: isequal 

        ! Auxiliary
        real(R8)                        :: tolfac

        ! Compute
        !========
        if (present(tolfacopt)) then 
            tolfac = tolfacopt*spacing(realval) 
        else 
            tolfac = spacing(realval) 
        end if 
        isequal = abs(realval-real(intval, kind=R8)) <= tolfac

    end function

end module