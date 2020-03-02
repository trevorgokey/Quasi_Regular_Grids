module grid_LJ_class

    use grid_base_class, only: grid_t, grid_base_init,initial_input_from_xyz

    type, extends( grid_t) :: grid_LJ_t
        integer :: opt_steps = 0, opt_tune = 0, opt_type = 0
        double precision :: opt_eps = 1d-7, cLJ = 1.0, mv_cutoff
    contains
        procedure, public :: V
        procedure, public :: Pair_LJ_NRG
        procedure, public :: U
        procedure, public :: optimize 
        !procedure :: P
    end type


contains

    function grid_LJ_init( input_fd, points_nr, grid_input, grid_from)
        implicit none
        type (grid_LJ_t), target :: grid_LJ_init
        !class (grid_t), intent(out), allocatable :: this
        integer,        intent(in ) :: input_fd, points_nr
        integer :: grid_input, opt_type, opt_steps,  opt_tune
        double precision :: opt_eps, cLJ, opt_dx
        character(len=*) :: grid_from
        !type (grid_LJ_t), target :: instance
        class (grid_t), pointer :: this
        namelist /grid_potential_LJ/ cLJ, opt_type, opt_steps, opt_eps, opt_tune, opt_dx
        read( unit=input_fd, nml=grid_potential_LJ)
        !allocate( grid_LJ_t::this)
        grid_LJ_init%cLJ = cLJ
        grid_LJ_init%opt_steps = opt_steps
        grid_LJ_init%opt_tune = opt_tune
        grid_LJ_init%opt_type = opt_type
        grid_LJ_init%opt_eps = opt_eps
        grid_LJ_init%mv_cutoff = opt_dx
        grid_LJ_init%points_nr = points_nr
        print*, "GRID Allocating", (3*points_nr + points_nr*points_nr) * 8, " bytes"
        !allocate( this%X(3, points_nr), this%Vij( points_nr, points_nr) )
        this => grid_LJ_init
        allocate( grid_LJ_init%X(3, points_nr), grid_LJ_init%Vij( points_nr, points_nr) )

        select case( grid_input)
            case (1)
            case (2)
                call this%initial_input_from_xyz( grid_from)
            case (3)
        end select
        !grid_LJ_init%V
        !this = grid_LJ_init
        !call grid_base_init( this, input_fd, points_nr, grid_input, grid_from)
        !grid_LJ_init%Vij = this%Vij
        !grid_LJ_init%X   = this%X
        !this = grid_LJ_init
        !grid_LJ_init = grid_LJ_init
        write (*,*) "GRID_LJ   init (grids alloc?)= ", allocated(grid_LJ_init%X) 
        write (*,*) "GRID_LJ       cLJ=", grid_LJ_init%cLJ, " points_nr ", grid_LJ_init%points_nr
    end function


    function V( this, x)
        class( grid_LJ_t), intent(in) :: this
        double precision :: V
        double precision, dimension(:), intent(in) :: x
        V = 0.0
        print*, "GRID_LJ   V" 
    end function

    function Pair_LJ_NRG(this, x1,x2)
        implicit none 
        class( grid_LJ_t) :: this
        double precision Pair_LJ_NRG, a, b, sigma1, sigma2
        double precision, dimension(:) :: x1, x2
        a=sum((x1(:)-x2(:))**2)
        Pair_LJ_NRG=0.0
        if ( a < 20**2 ) then
        sigma1=this%cLJ * (this%system%P( x1) * this%points_nr)**(-1./this%system%d)    
        sigma2=this%cLJ * (this%system%P( x2) * this%points_nr)**(-1./this%system%d)    
        b=(sigma2**2/a)**3
        a=(sigma1**2/a)**3
        Pair_LJ_NRG=a**2-a+b**2-b
        !print*,"GRID_LJ::NRG sigma1,2=", sigma1, sigma2, " a,b ", a, b, " E=", Pair_LJ_NRG
        end if
    end function Pair_LJ_NRG

    function U( this)

        class( grid_LJ_t) :: this
        integer :: i, j, points_nr
        double precision, dimension(this%points_nr,this%points_nr)  :: Vij
        double precision, dimension(3,this%points_nr)  :: X
        Vij = this%Vij
        X   = this%X
        points_nr = this%points_nr

        !$OMP parallel default(none) shared(this,points_nr) private(i,j)
        !$OMP do 
        !print*, "GRID_LJ::U points_nr", points_nr
        do i=2,this%points_nr
            do j=1,i-1
                !write(*,*) "GRID_LJ::U ", i,j,Npoints,i-1, this%X(:,i), this%X(:,j)
                this%Vij(i,j)=this%Pair_LJ_NRG( this%X(:,i), this%X(:,j))
                this%Vij(j,i)=this%Vij(i,j)
            enddo
        enddo
        !$OMP end do
        !$OMP end parallel
    end function

    subroutine optimize( this)
        class (grid_LJ_t) :: this
        type (grid_LJ_t)  :: obj

        double precision, dimension(:), allocatable :: s, x0
        double precision, dimension(:), allocatable :: U_move
        double precision :: Delta_E, mv_cutoff=1.0, deltae1=0.0, V, P
        integer :: i, j, k, counter=0, Npoints, accept=0, plt_count=0, MMC_freq, N_MC, newmove=0
        !double precision, dimension(this%points_nr,this%points_nr)  :: Vij
        !double precision, dimension(3,this%points_nr)  :: X


        !Vij = this%Vij
        !X   = this%X
        
        obj = this

        print*, "GRID_LJ OPTIMIZER", this%opt_tune, this%cLJ

        MMC_freq = this%opt_tune
        N_MC     = this%opt_steps
        mv_cutoff= this%mv_cutoff

        !open(unit=18,file='mv_cut.dat')
        !open(unit=19,file='delE.dat')

        Npoints = this%points_nr

        allocate( U_move(this%points_nr), s(this%system%d), x0(this%system%d))

        do i=1,N_MC
        newmove = 0
        write(*,*) i,N_MC
        do k=1, Npoints
        !==============================================================================!
        !                       Randomly Select Atom to Move
        !==============================================================================!
            call random_number(s) 
            !k=random_integer(1,Npoints, s(1))
        !==============================================================================!
        !                           Generate trial move
        !        random numbers generated (0,1), make it (-1,1) ==> s=2*s-1
        !==============================================================================!
            x0=this%X(:,k)+mv_cutoff*(2*s-1)
        !==============================================================================!
        !                   Only consider point if V(trial) < Ecut 
        !                   Compute Energy Change due to Trial Move
        !==============================================================================!
            !print*, "x0 is", x0
            V = this%system%V( x0)
            if( V < this%system%E_cut) then
                counter=counter+1
                Delta_E=0d0
                P = this%system%P(x0)
                U_move(k)= P
                !print*, "Umove prior to loop is ", U_move(k), "allocated?", allocated(U_move), " Valloc? ", allocated(this%Vij)
        !  !$OMP parallel default(none) shared(Npoints,this%Vij,this%X,k,U_move,this%system%X,Natoms,x0) private(j,diff) reduction(+:Delta_E)
        !  $OMP p#arallel default(none) shared(Npoints,Vij,X,k,U_move,x0) private(j,diff) reduction(+:Delta_E)
        !$OMP parallel default(none) shared(Npoints,this,k,U_move,x0) private(j,diff) reduction(+:Delta_E)
        !$OMP do 
                do j=1,Npoints
                    !write(*,*) "k in parallel is", k, "j is", j
                    if(j.ne.k) then
                        U_move(j)=this%Pair_LJ_NRG( this%X(:,j), x0)
                        Delta_E = Delta_E+ ( this%Vij(j,k) - U_move(j))
                        !write(*,*) j,k,"Ujk=", this%Vij(j,k), " DE= ", Delta_E, " trial U_move was", U_move(j), &
                        !    & "leading to ", this%Vij(j,k)-U_move(j)
                    endif
                enddo
        !$OMP end do
        !$OMP end parallel
                write(*,*) "Delta_E from reduce is", Delta_E
        !==============================================================================!
        !                  Accept any trial that decreases the energy 
        !==============================================================================!
                write(*,*) "new is", U_move(k), "old was", this%Vij(k,k), "P(x) is", P
                if(Delta_E>0d0)then
                    this%Vij(:,k)=U_move(:)
                    this%Vij(k,:)=U_move(:)
                    accept=accept+1
                    newmove = 1
                    this%X(:,k)=x0(:)
                    deltae1=deltae1+Delta_E
                    write(*,*) k, Delta_E, "Accepted! Displacement was", &
                        & sqrt(sum((mv_cutoff*(2*s-1))**2)), "V(x0)", V


                endif
            endif
        enddo
        if( newmove > 0) then
            open(unit=20,file='grid.xyz', access='APPEND')
        
                write(20,*) Npoints
                write(20,*)
                do j=1,Npoints
                    write(20,*) "C     ", this%X(:,j)
                enddo
            close(20)
        end if
    !==============================================================================!
    !                           Update Cutoff Paramater
    !           acceptance rate ~50%, adjust trial displacement accordingly
    !==============================================================================!
        if(mod(i,MMC_freq)==0)then
            write(*,*) 'MMC Iteration', i
            if(dble(accept)/counter<0.5)then 
                mv_cutoff=mv_cutoff*0.9
            else
                mv_cutoff=mv_cutoff*1.1
            endif
            if (mv_cutoff < .01) then
                mv_cutoff = .01
            endif
            accept=0
            counter=0
            !call Moments(Moment,x)
            !write(*,*) 'Moments:'
            !write(*,*) Moment(1:9)
            write(*,*) 'mv cutoff', mv_cutoff
            write(*,*) 'Deltae1==>', deltae1
            plt_count=plt_count+1
            write(18,*) plt_count, mv_cutoff
            write(19,*) plt_count, deltae1
            deltae1=0
        endif
        !write(*,*) "wrote xyz"
        !open(unit=21,file='grid.dat')
        !    do j=1,Npoints
        !        write(21,*) x(:,j)
        !    enddo
        !close(21)
    enddo
    close(18)
    close(19)

    end subroutine

end module grid_LJ_class
