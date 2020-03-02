
module grid_base_class
    use space_base_class
    use system_base_class, only: system_t
    implicit none

    type, extends(space_t), abstract :: grid_t
        double precision, allocatable, dimension(:,:) :: X
        double precision, allocatable, dimension(:,:) :: Vij
        integer :: points_nr = 0
        class(system_t), allocatable :: system
    contains
        procedure :: grid_set_system
        procedure :: initial_input_from_xyz
        procedure(optimize), deferred, public :: optimize
    end type 

    abstract interface
        subroutine optimize( this)
            import :: grid_t
            implicit none
            class (grid_t) :: this
        end subroutine
    end interface

contains

    subroutine grid_base_init( this, input_fd, N, grid_input, grid_from)
        class (grid_t), intent(out), pointer :: this
        integer,        intent(in ) :: input_fd, N
        integer :: grid_input
        character(len=*) :: grid_from
        call space_base_init( this, input_fd)
        !this%points_nr = N
        !allocate( this%X(3, N), this%Vij( N, N) )
        print*, "GRID init"
        select case( grid_input)
            case (1)
            case (2)
                call this%initial_input_from_xyz( grid_from)
            case (3)
        end select
        !this%test = input_fd
    end subroutine


    subroutine grid_set_system( this, system)
        class (grid_t)   :: this
        class (system_t) :: system
        this%system = system
        write (*,*) "GRID grid_set_system"
    end subroutine

    subroutine initial_input_from_xyz(this, fname )
        integer :: i, j, k, readin, fid=15, Npoints, Natoms
        class (grid_t) :: this
        character(len=*) :: fname
        Npoints = this%points_nr
        if( allocated(this%X) ) then
            print*, "ALLOCATED"
        else
            print*, "NOT ALLOCATED...."
        end if
        open( fid, file=fname)
        read( fid, *) Natoms
        print*, "GRID Reading initial crds from ", trim(fname), " Natoms= ", Natoms, "Npoints=", Npoints
        if(Natoms < Npoints) then
            print*, "This file has fewer points than requested"
            stop
        end if
        read( fid, *)
        k = 0
        do i = 1, Npoints
        k = k + 1
        if (k <= Natoms) then
            read(fid, *) this%X(:,i)
            readin = readin + 1
            write(*,*) "Read i=", i, "k=",k,Natoms, Npoints,"ik=",(i-1)*Natoms/Npoints, this%X(:,i)
            do j=1 , Natoms/Npoints - 1
                k = k+1
                if (k < Natoms) then
                    !write(*,*) "skipping", i, j, k, "steps", nint(dble(Natoms)/dble(Npoints)) - 1
                    read(fid, *)
                else if ( k == Natoms) then
                    write(*,*) "keeping last", i, j, k, "steps", nint(dble(Natoms)/dble(Npoints)) - 1
                    read(fid, *) this%X(:,i+1)
                    readin = readin + 1
                else
                    write(*,*) "skipping", i, j, k
                endif 
            enddo
        endif
        end do
        close( fid)
    end subroutine

end module grid_base_class
