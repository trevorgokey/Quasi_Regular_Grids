module system_basis_class

    use system_base_class, only: system_t, system_base_init
    implicit none

    type, extends( system_t), abstract :: system_basis_t
        character(len=512) :: centers_fnm
        double precision, allocatable, dimension(:,:) :: centers
        integer :: Ncenters
    contains
        private
        procedure, public  :: read_centers
    end type

contains

    subroutine system_basis_init( this, fid)
        class (system_t), intent(out) :: this
        integer,          intent(in ) :: fid
        call system_base_init( this, fid)
        print*, "SYSTEM_BASIS init"
    end subroutine

    subroutine read_centers( this)
        implicit none
        class (system_basis_t) :: this
        integer :: i, fid=15
        this%d = 3
        open( fid, file=this%centers_fnm, iostat=i, action='read')
        if (i .ne. 0) then
            print *, "Could not open centers file: ", trim(this%centers_fnm)
            stop
        end if
        read( fid, *) this%Ncenters
        allocate( this%X( this%d, this%Ncenters))
        read( fid, *)
        do i = 1, this%Ncenters
            read( fid, *) this%X(:, i)
        end do
        close( fid)
        print*, "SYSTEM_BASIS read_centers file ", trim(this%centers_fnm), " Ncenters=", this%Ncenters, "d=", this%d 
        print*, "coords"
        print*, this%X
    end subroutine

end module


