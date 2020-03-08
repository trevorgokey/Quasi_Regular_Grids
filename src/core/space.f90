
module space_base_class

    implicit none

    type, abstract :: space_t
        integer :: d
    contains
        private
        procedure(V), deferred, public :: V
        procedure(dVdx), deferred, public :: dVdx
    end type 

    abstract interface
        function V( this, x)
            import :: space_t
            implicit none
            class (space_t), intent(in) :: this
            double precision :: V
            double precision, dimension(:), intent(in) :: x
        end function

        function dVdx( this, x)
            import :: space_t
            implicit none
            class (space_t), intent(in) :: this
            double precision :: V
            double precision, dimension(this%d), intent(in) :: x
            double precision, dimension(this%d) :: dVdx
        end function
    end interface



contains

    subroutine space_base_init( this, input_fd)
        class (space_t), intent(out) :: this
        integer,         intent(in ) :: input_fd
        print*, "SPACE_BASE init"
    end subroutine

end module space_base_class
