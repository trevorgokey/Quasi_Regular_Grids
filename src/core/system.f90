
module system_base_class

    use space_base_class, only: space_t, space_base_init
    implicit none
    type, extends(space_t), abstract :: system_t
        double precision, allocatable, dimension(:,:) :: X
        double precision :: E_cut
        integer :: d=3
    contains
        procedure, public :: P
    end type 

contains

    subroutine system_base_init( this, fid)
        class (system_t), intent(out) :: this
        integer,          intent(in ) :: fid
        print*, "SYSTEM init"
    end subroutine

    function P( this, x)
        implicit none 
        class (system_t) :: this
        double precision, dimension(:) :: x
        double precision::P,integral_P=1.0

        P = this%V( x)
        if( P < this%E_cut) then
           P = ( this%E_cut - P)/integral_P          
        else                                             !set equal to 0 if beyond Ecut
           P = 1d-20
        end if
    end function P

end module system_base_class
