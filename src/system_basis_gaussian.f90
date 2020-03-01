module system_basis_gaussian_class

    use system_base_class,  only: system_t
    use system_basis_class, only: system_basis_t, system_basis_init
    implicit none

    type, extends( system_basis_t) :: system_basis_gaussian_t
        double precision :: cPrefac, cSigma

    contains
        private
        procedure, public :: V
    end type

contains

    subroutine system_basis_gaussian_init( this, fid)
        class (system_t), allocatable, intent(out) :: this
        integer,                       intent(in ) :: fid
        type (system_basis_gaussian_t) :: instance
        double precision :: cPrefac, cSigma, E_cut
        character( len=512) :: centers
        integer :: post_proc

        namelist /system_potential_gaussian_basis/ cPrefac, cSigma, centers, post_proc, E_cut
        read( unit=fid, nml=system_potential_gaussian_basis)
        write (*,*) "SYSTEM_BASIS_GAUSS   init"
        allocate( system_basis_gaussian_t::this)
        instance%cPrefac     = cPrefac
        instance%cSigma      = cSigma
        instance%centers_fnm = trim(centers)
        call instance%read_centers()

        call system_basis_init( this, fid)
        this = instance
        this%E_cut = E_cut

    end subroutine


    function V( this, x)
        implicit none
        class (system_basis_gaussian_t), intent(in) :: this
        integer::i
        double precision :: V, dist
        double precision, dimension(:), intent(in) :: x
        
        V = 0.0
        !print*, "SYS_GBAS   Ncenters", this%Ncenters
        do i = 1, this%Ncenters
            dist = (sum( (x(:) - this%X(:,i))**2))
            !if ( dist < sigma*3 ) then
            V = V - this%cPrefac*exp(-dist/this%cSigma) 
            !print*, "SYS_GBAS   energy V dist x X", V, dist, x(:), this%X(:,i)
            !print*, "SYS_GBAS   energy dist", sqrt(sum( (x(:) - this%X(:,i))**2)) !x(:) - this%X(:,i)
            !end if
        end do
        !print*, "SYS_GBAS   energy", V, dist, x(:)! , this%X(:,i)
    end function



end module
