
module system_base_class

    use space_base_class, only: space_t, space_base_init
    implicit none
    type, extends(space_t), abstract :: system_t
        double precision, allocatable, dimension(:,:) :: X
        double precision :: E_cut
    contains
        procedure, public :: P
        procedure, public :: dPdx
    end type 

contains

    subroutine system_base_init( this, fid)
        class (system_t), intent(out) :: this
        integer,          intent(in ) :: fid
        this%d = 3
        print*, "SYSTEM init"
    end subroutine

    function P( this, x)
        implicit none 
        class (system_t) :: this
        integer :: d
        double precision, dimension(:) :: x
        double precision::P,integral_P=1.0

        P = this%V( x)
        !print*,"P:: V is ", P, "Ecut - V= ", this%E_cut - P
        d = this%d
        if( this%E_cut - P > 0.0 ) then
           P = ( this%E_cut - P)**(d/2.0) 
        else
           P = exp(this%E_cut - P)
        endif
    end function P

    function dPdx( this, x)
        implicit none 
        class (system_t) :: this
        double precision, dimension(this%d) :: x
        double precision::integral_P=1.0, V, d, diff
        double precision, dimension(this%d) :: dPdx

        V = this%V( x)
        d = this%d
        diff = this%E_cut - V
        if( diff < 0.0) then
            dPdx(:) = -this%dVdx( x)
        else
            !print*, "dPdx was ", dPdx
            dPdx(:) = -(d/2.0)*( diff)**(d/2.0 - 1) * this%dVdx( x)
            !print*, "E-V= ", ( this%E_cut - V), " power= ", ( this%E_cut - V)**(d/2.0 - 1), " dVdx ",  this%dVdx( x)
            !dPdx(:) = dPdx(:) * exp(this%E_cut - V)*this%dVdx( x)
            !print*, "dPdx is now", dPdx
        endif
    end function dPdx

end module system_base_class
