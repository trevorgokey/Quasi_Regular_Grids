module grid_LJ_class

    use grid_base_class, only: grid_t, grid_base_init,initial_input_from_xyz

    type, extends( grid_t) :: grid_LJ_t
        integer :: opt_steps = 0, opt_tune = 0, opt_type = 0, saveevery = 1
        double precision :: opt_eps = 1d-7, cLJ = 1.0, mv_cutoff
        logical :: adaptive  = .false.
    contains
        private
        procedure, public :: V
        procedure, public :: Pair_LJ_NRG
        procedure, public :: U
        procedure, public :: optimize 
        procedure, public :: dVdx
        procedure :: duijdx
        procedure, public :: d2Vdx2
        !procedure :: P
    end type


contains

    function grid_LJ_init( input_fd, points_nr, grid_input, grid_from)
        implicit none
        type (grid_LJ_t), target :: grid_LJ_init
        !class (grid_t), intent(out), allocatable :: this
        integer,        intent(in ) :: input_fd, points_nr
        integer :: grid_input, opt_type, opt_steps,  opt_tune, saveevery
        double precision :: opt_eps, cLJ, opt_dx
        character(len=*) :: grid_from
        !type (grid_LJ_t), target :: instance
        class (grid_t), pointer :: this
        namelist /grid_potential_LJ/ cLJ, opt_type, opt_steps, opt_eps, opt_tune, opt_dx, saveevery
        read( unit=input_fd, nml=grid_potential_LJ)
        !allocate( grid_LJ_t::this)
        grid_LJ_init%cLJ = cLJ
        grid_LJ_init%opt_steps = opt_steps
        grid_LJ_init%opt_tune = opt_tune
        grid_LJ_init%opt_type = opt_type
        grid_LJ_init%opt_eps = opt_eps
        grid_LJ_init%mv_cutoff = opt_dx
        grid_LJ_init%points_nr = points_nr
        grid_LJ_init%saveevery = saveevery
        print*, "GRID Allocating", (3*points_nr + points_nr*points_nr) * 8, " bytes"
        !allocate( this%X(3, points_nr), this%Vij( points_nr, points_nr) )
        this => grid_LJ_init
        if( .not. allocated(grid_LJ_init%X) ) then
            allocate( grid_LJ_init%X( 3, this%points_nr))
        endif
        if( .not. allocated(grid_LJ_init%Vij) ) then
            allocate( grid_LJ_init%Vij( points_nr, points_nr) )
        endif

        select case( grid_input)
            case (1)
                call this%initial_input_from_xyz( grid_from, 1)
            case (2)
                call this%initial_input_from_xyz( grid_from, 0)
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
        double precision Pair_LJ_NRG, a, b, sigma1, sigma2, P
        double precision, dimension(:) :: x1, x2
        a=sum((x1(:)-x2(:))**2)
        Pair_LJ_NRG=0.0
        if( a > 0.0) then
            !if ( a < 20**2 ) then
            P = this%system%P( x1)
            !print*,"GRID_LJ::NRG P(x1)= ", P
            sigma1=this%cLJ * (P * this%points_nr)**(-1./this%system%d)
            if(P == 0.0) sigma1 = 0.0
            P = this%system%P( x2)
            !print*,"GRID_LJ::NRG P(x2)= ", P
            sigma2=this%cLJ * (this%system%P( x2) * this%points_nr)**(-1./this%system%d)    
            if(P == 0.0) sigma2 = 0.0
            b=(sigma2**2/a)**3
            a=(sigma1**2/a)**3
            Pair_LJ_NRG=a**2-a+b**2-b
            !end if
        endif
        !print*,"GRID_LJ::NRG sigma1,2=", sigma1, sigma2, " a,b ", a, b, " E=", Pair_LJ_NRG, " d= ", sum((x1(:)-x2(:))**2)
        !print*,"GRID_LJ::NRG x1 ", x1, " x2 ", x2
        !print*,"GRID_LJ::NRG P1 ", this%system%P( x1), " P2 ", this%system%P( x2)
    end function Pair_LJ_NRG

    function U( this)

        class( grid_LJ_t) :: this
        integer :: i, j, points_nr
        double precision, dimension(this%points_nr,this%points_nr)  :: Vij
        double precision, dimension(3,this%points_nr)  :: X
        double precision :: U
        Vij = this%Vij
        X   = this%X
        points_nr = this%points_nr

        !$OMP parallel default(none) shared(this) private(i,j,points_nr)
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
        U = sum(this%Vij)/2.0
    end function

    function duijdx(this, x1, x2)
        implicit none 
        class( grid_LJ_t) :: this
        integer :: i, j, N, d, k
        double precision r, sigma1, dx, a, b, P, m
        double precision, dimension(this%system%d) :: x1, x2, duijdx, dP, dV
        duijdx(:) = 0.0
        r=sqrt(sum((x1(:)-x2(:))**2))
        !if ( a < 20**2 ) then
        N = this%points_nr
        d = this%system%d


        !print*, "x1 x2", x1, x2
        P = this%system%P( x1)
        dP = this%system%dPdx( x1)
        m = 12
        k = 6
        do i=1, size(x1)
            dx = x2(i) - x1(i)
            !print*, "dx ", dx
            !print*, "dPdx for i=", i, P, dP
            a = 0.0
            b = 0.0
            !print*, "i j sigma1", i, j, sigma1
            !print*,"a=", a, " sigma1 ", sigma1, " P ", P

            ! uji case
            if( P .ne. 0.0) then
                sigma1 = this%cLJ * (P * N)**(-1./d) 
                !a = dx/r**2 - 1.0/(P*d) * dP(i)
                a = dx/r**2 - 1.0/(P*d) * dP(i)
                b = m*(sigma1/r)**m - k*(sigma1/r)**k
                duijdx(i) = duijdx(i) + a*b
            endif
        enddo
        
        !print*,"GRID_LJ::NRG sigma1,2=", sigma1, sigma2, " a,b ", a, b, " E=", Pair_LJ_NRG
    end function duijdx

    function dVdx( this, x)

        class( grid_LJ_t), intent(in) :: this
        double precision, dimension(this%system%d), intent(in)  :: x
        integer :: i, j, k, points_nr
        double precision, dimension(this%points_nr,this%points_nr)  :: Vij
        double precision, dimension(this%system%d)  :: dVdx
        double precision :: r = 0.0

        !print*, "GRID_LJ::U points_nr", points_nr
        dVdx(:) = 0.0
        do i=1,this%points_nr
            !print*,"  j", i, " dist=", sqrt(sum((this%X(:,i) - x)**2))
            if( sqrt(sum((this%X(:,i) - x)**2)) > 0.0) then
                ! this version pushes everything out ( center less dense
                ! duij did -dx in dudxj term, sum both terms
                ! + x X - X x
                !dVdx(:) = dVdx(:)   - this%duijdx( x, this%X(:,i)) + this%duijdx( this%X(:,i), x) !

                ! this clumps
                !dVdx(:) = dVdx(:)   + this%duijdx( x, this%X(:,i)) - this%duijdx( this%X(:,i), x) !

    
                ! -  repels everything outside cut, + brings everything to origin
                !dVdx(:) = dVdx(:)   + this%duijdx( x, this%X(:,i)) !- this%duijdx( this%X(:,i), x) !

                ! -  repels, but everything seems more coupled 
                !dVdx(:) = dVdx(:)   - this%duijdx( this%X(:,i), x) !

                ! this + - produces the bee effect even when P1 and P2
                !dVdx(:) = dVdx(:)   + this%duijdx( this%X(:,i), x) !+ this%duijdx( x, this%X(:,i))

                ! + does the repulsive bee thing, - attracts
                !dVdx(:) = dVdx(:)   + this%duijdx( this%X(:,i), x) !+ this%duijdx( x, this%X(:,i))

                ! - repels each other, + does a local clumping
                !dVdx(:) = dVdx(:)   - this%duijdx( x, this%X(:,i))
                !print*,"duijdx i j", x, i, this%duijdx( x, this%X(:,i)), "total=", dVdx(:) 
                dVdx(:) = dVdx(:) + this%duijdx( x, this%X(:,i)) !+ this%duijdx( x, this%X(:,i))
            else
                !dVdx(:) = dVdx(:) + this%system%dVdx( x)+ this%system%dVdx( this%X(:,i))
                ! + is repulsive
                dVdx(:) = dVdx(:) - ((this%system%dPdx( x)))! + this%system%dPdx( this%X(:,i)))
            endif
            !print*,"  j", i, " dist=", sqrt(sum((this%X(:,i) - x)**2)), " dUdx", dVdx

                !write(*,*) "GRID_LJ::U ", i,j,Npoints,i-1, this%X(:,i), this%X(:,j)
        enddo
    end function

    subroutine optimize( this)
        class (grid_LJ_t) :: this
        real, dimension(this%system%d,this%points_nr) :: buf
        integer :: steps
        open(unit=20,file='grid.xyz', status='REPLACE')
        close(20)
        this%adaptive=.true.
        call optimize_MCCG( this)
        call optimize_CG( this)
        call optimize_MCCG( this)
    end subroutine

    subroutine optimize_CG( this)
        class (grid_LJ_t) :: this

        !double precision, dimension(:), allocatable :: s, x0
        !double precision, dimension(:), allocatable :: U_move
        double precision, dimension(this%system%d,this%points_nr)  :: dUdx, totaldUdx, sdVdx
        double precision, dimension(this%points_nr)  :: buff
        double precision :: Delta_E, mv_cutoff=1.0, mv_cutoff_ref=1.0, deltae1=0.0, Eprev, Vprev, V, V_trial, P, E, SDX, dx=1.0
        double precision :: Vs, Vsprev, F, Fprev, time1, time2
        integer :: i, j, k, counter=0, Npoints, accept=0, hotaccept=0,plt_count=0, MMC_freq, N_MC, newmove=0, framesaved=1
        integer :: saveevery
        logical :: converged = .false., adaptive = .false.
        character(1) :: saved = 'N'
        logical, parameter :: timing = .false.
        !double precision, dimension(this%points_nr,this%points_nr)  :: Vij
        !double precision, dimension(3,this%points_nr)  :: X
        !character(len=*), parameter :: fmt = "(I8,A1,I10,A5,F6.3,A5,E8.6,A5,E8.6,A5,E8.6,5A,E8.6,A7,F6.3,5A,I4,I6)"
        character(len=*), parameter :: fmt = '(I10, " /", I10,1x, I10, "  F=", ES15.8, "  DF=", ES15.8, "  E=", ES15.8, &
            "   DE=", ES15.8, "   DX=", ES15.8, "   V=", ES15.8, "   dV=", ES15.8,"   Vg=", ES15.8, "   dVg= ", ES15.8 )'

        !Vij = this%Vij
        !X   = this%X
        

        print*, "GRID_LJ OPTIMIZER", this%opt_tune, this%cLJ, " steps=", N_MC

        MMC_freq = this%opt_tune
        N_MC     = this%opt_steps
        mv_cutoff= this%mv_cutoff
        mv_cutoff_ref= this%mv_cutoff
        saveevery = this%saveevery 
        dx = .01 !this%mv_cutoff
        adaptive = this%adaptive

        !open(unit=18,file='mv_cut.dat')
        !open(unit=19,file='delE.dat')

        Npoints = this%points_nr

        !allocate( U_move(this%points_nr), s(this%system%d), x0(this%system%d))


        i = 0
        open(unit=20,file='grid.xyz', access='APPEND')
            write(20,*) Npoints
            write(20,*) "step", i, "frame", framesaved
            do j=1,Npoints
                write(20,*) "C     ", this%X(:,j)
            enddo
        close(20)
        !print*,"X allocated?", allocated( this%X)
        Eprev = 0.0
        Vprev = 0.0
        Vgprev = this%U() 
        totaldUdx(:,:) = 0.0
        sdVdx = 0.0
        buff(:) = 0.0
        if( timing) print*, "SYSTEM V"
        call cpu_time( time1)
        !$OMP parallel default(none) shared(this, buff, Npoints) private(k)
        !$OMP do 
        do k=1,Npoints
            buff(k) = this%system%V( this%X(:,k))
            !do j=k+1, Npoints
            !    !print*, "X k j=", j, k, this%X(:,k), this%X(:,j)
            !    Vprev = Vprev + this%Pair_LJ_NRG( this%X(:,j), this%X(:,k))
            !    !print*, "V j k=", j, k, Vprev
            !enddo
        enddo
        !$OMP end do
        !$OMP end parallel

        do k=1,Npoints
            Vprev = Vprev + buff(k)
        enddo
        call cpu_time( time2)
        if( timing) print*, "TIME ", time2 - time1

        !print*,"grid ", this%X
        do
            i = i + 1
            if ( N_MC == 0 .or. ( N_MC > -1 .and. i == N_MC)) then
                exit
            endif
            !print*,"X"
            !print*,this%X
            !do k=1,Npoints
            !    do j=1, k-1
            !        !print*, "X k j=", j, k, this%X(:,k), this%X(:,j)
            !        this%Vij(j,k) = this%Pair_LJ_NRG( this%X(:,j), this%X(:,k))
            !        !print*, "V j k=", j, k, this%Vij(j,k)
            !    enddo
            !    do j=1, k-1
            !        this%Vij(k,j) = this%Vij(j,k)
            !        !print*, "V k j=", k, j, this%Vij(k,j)
            !    enddo
            !    this%Vij(k,k) = 0.0
            !enddo
            !print*, "V="
            !print*, this%Vij
            newmove = 0
            !print*, "eval grad", i
            F = 0.0
            if( timing) print*,"DVDX"
            call cpu_time( time1)

            if ( .true.) then
                !$OMP parallel default(none) shared(this,dUdx,sdVdx,buff, mv_cutoff, Npoints) private(j,E, F)
                !$OMP do 
                do j=1, Npoints
                    dUdx(:,j) = -this%dVdx( this%X(:, j))
                    !print*,j, "dUdx is", dUdx(:,j)

                    F = sqrt(sum(dUdx(:,j)**2))
                    buff(j) = F
                    if( F > mv_cutoff ) then
                        dUdx(:,j) = dUdx(:,j) / F * mv_cutoff
                    endif
                enddo
                !$OMP end do
                !$OMP end parallel
                F = 0.0
                do j=1, Npoints
                    F = F + buff(j)
                enddo
                F = F/Npoints
            else
                !$OMP parallel default(none) shared(this,dUdx,sdVdx,buff, mv_cutoff, Npoints) private(j,E, F)
                !$OMP do 
                do j=1, Npoints
                    dUdx(:,j) = -this%dVdx( this%X(:, j))
                    sdVdx(:,j) = -this%system%dPdx( this%X(:, j))
                enddo
                !$OMP end do
                !$OMP end parallel
                !print*, "distances ", sqrt(sum(dUdx**2)), sqrt(sum(sdVdx**2))
                dUdx = dUdx + sdVdx 
                dUdx = dUdx / 2.0
                !dx = sum(dUdx/sqrt(sum(dUdx**2)) * sdVdx/sqrt(sum(sdVdx**2)))
                !if (sqrt(sum(dUdx**2)) > sqrt(sum(sdVdx**2))) then
                !    dUdx(:,:) = sdVdx
                !    dx = dx * sqrt(sum(sdVdx**2))
                !else
                !    dx = dx * sqrt(sum(dUdx**2))
                !endif
                !dx=dx * min( sqrt(sum(dUdx**2)), )
                !if (dx > 0) then

                !    dUdx = (dUdx + sdVdx) / 2.0
                !    dx = sqrt(sum(dUdx**2))
                !    print*,"QUIT! dx= ", dx, "F=", F
                !end if
                F = sqrt(sum(dUdx**2))
                dx=min(dx, mv_cutoff)
                dUdx = dUdx / sqrt(sum(dUdx**2)) * dx
            endif



            call cpu_time( time2)
            if( timing) print*, "TIME ", time2 - time1
            if( timing) print*,"SUM DVDX"
            call cpu_time( time1)
            !$OMP parallel default(none) shared(dUdx, buff, Npoints) private(j)
            !$OMP do 
            do j=1, Npoints
                buff(j) = sqrt(sum(dUdx(:,j)**2))
            enddo
            !$OMP end do
            !$OMP end parallel

            E = 0.0
            do j=1, Npoints
                E = E + buff(j)
            enddo
            E = E/Npoints
            call cpu_time( time2)
            if( timing) print*, "TIME ", time2 - time1

            if( timing) print*,"DVDX *DX"
            call cpu_time( time1)
            !$OMP parallel default(none) shared(dUdx, this, Npoints) private(j)
            !$OMP do 
            do j=1, Npoints
                this%X(:,j) = this%X(:,j) + dUdx(:,j)
            enddo
            !$OMP end do
            !$OMP end parallel
            call cpu_time( time2)
            if( timing) print*, "TIME ", time2 - time1
            
            V = 0.0
            if( timing) print*,"U"
            call cpu_time( time1)
            Vg = this%U()
            call cpu_time( time2)
            if( timing) print*, "TIME ", time2 - time1

            if( timing) print*,"SYSTEM V"
            call cpu_time( time1)
            !$OMP parallel default(none) shared(this, buff, Npoints) private(j)
            !$OMP do 
            do j=1,Npoints
                buff(j) = this%system%V( this%X(:,j)) 
                !do j=k+1, Npoints
                !    !print*, "X k j=", j, k, this%X(:,k), this%X(:,j)
                !    V = V + this%Pair_LJ_NRG( this%X(:,j), this%X(:,k))
                !    !print*, "V j k=", j, k, this%Vij(j,k)
                !enddo
            enddo
            !$OMP end do
            !$OMP end parallel

            V = 0.0
            do j=1,Npoints
                V = V + buff(j) 
                !do j=k+1, Npoints
                !    !print*, "X k j=", j, k, this%X(:,k), this%X(:,j)
                !    V = V + this%Pair_LJ_NRG( this%X(:,j), this%X(:,k))
                !    !print*, "V j k=", j, k, this%Vij(j,k)
                !enddo
            enddo
            !newmove=1
            !print fmt,  i, N_MC, 'N', F, F - Fprev, V+Vg, V+Vg-Vprev-Vgprev, E, V, V - Vprev, Vg, Vg - Vgprev
            call cpu_time( time2)
            if( timing) print*, "TIME ", time2 - time1
            if( adaptive) then
            if( i > 1) then
                !if( sign(1.0d1,Eprev) .ne. sign(1.0d1,E) .or. V-Vprev > 0) then

                if( V+Vg-Vprev-Vgprev > 0 ) then
                    if( mv_cutoff .eq. mv_cutoff_ref) then
                        totaldUdx(:,:) = this%X(:,:) 
                    endif
                        
                    mv_cutoff = mv_cutoff * .95
                    newmove=0
                    !$OMP parallel default(none) shared(this, dUdx, Npoints) private(j)
                    !$OMP do 
                    do j=1, Npoints
                        this%X(:,j) = this%X(:,j) - dUdx(:,j)
                    enddo
                    !$OMP end do
                    !$OMP end parallel
                else
                    !print*, "step", i, "USUM", sum( dUdx), "DX", E/Npoints, " V= ", V, " dV= ", V - Vprev
                    !print fmt,  i, N_MC, '!', F, F - Fprev, V+Vg, V+Vg-Vprev-Vgprev, E, V, V - Vprev, Vg, Vg - Vgprev
                    newmove=1
                    !Eprev = E
                    mv_cutoff = mv_cutoff_ref
                    !Vprev = V
                    !Vgprev = Vg
                    !Fprev = F
                endif
                if (newmove .eq. 0 .and. mv_cutoff < this%opt_eps*10 .and. mv_cutoff_ref > this%opt_eps) then
                    this%X(:,:) = totaldUdx(:,:) 
                    newmove=1
                    mv_cutoff_ref = mv_cutoff_ref * .95
                    !print fmt,  i, N_MC, '!', F, F - Fprev, V+Vg, V+Vg-Vprev-Vgprev, E, V, V - Vprev, Vg, Vg - Vgprev
                    !print*, "**** mv_cutoff shrunk to", mv_cutoff_ref
                    mv_cutoff = mv_cutoff_ref
                    !Eprev = E
                    !Vprev = V
                    !Vgprev = Vg
                    !Fprev = F
                    continue
                endif
            endif
            else
                newmove=1
            endif
            !if( E/Npoints < this%opt_eps .or. mv_cutoff_ref < this%opt_eps) then
            if( newmove > 0) then
                converged = .false.
                if( abs(F-Fprev) < this%opt_eps) then
                    print *, "force convergence" 
                    converged = .true.
                endif
                if( abs(V+Vg - Vprev-Vgprev) < this%opt_eps) then
                    print *, "energy convergence"
                    converged = .true.
                endif
                if( mod( i, saveevery) == 0 .or. converged) then
                    print fmt,  i, N_MC, framesaved, F, F - Fprev, V+Vg, V+Vg-Vprev-Vgprev, E, V, V - Vprev, Vg, Vg - Vgprev
                    !print*, "step", i, " F= ", F, " dF= ", F - Fprev, " E= ", V+Vg, "DX", E, " V= ", V, " dV= ", V - Vprev, &
                    !    & " Vg= ", Vg, " dVg= ", Vg - Vgprev
                    open(unit=20,file='grid.xyz', access='APPEND')
                        write(20,*) Npoints
                        write(20,*) "step", i, "frame", framesaved, " E ", E, " V ", V, " Vg ", Vg, " F ", F
                        do j=1,Npoints
                            write(20,*) "C     ", this%X(:,j)
                        enddo
                    close(20)
                    saved='Y'
                    newmove = 1
                    framesaved = framesaved + 1
                endif
                if( converged) then
                    exit
                endif
                Eprev = E
                Vprev = V
                Vgprev = Vg
                Fprev = F
            endif

!            newmove = 0
!            !write(*,*) i,N_MC
!            E = sum(this%Vij) / 2
!            V = 0.0
!            SDX = 0.0
!            do k=1, Npoints
!                V = V + this%system%V( this%X(:,k))
!            enddo
!            do k=1, Npoints
!
!                call random_number(s) 
!                x0=this%X(:,k)+mv_cutoff*(2.0*s-1.0)
!                !print*, "DELTA ", mv_cutoff*(2.0*s-1.0)
!                !print*, "x0 is", x0
!                V_trial = this%system%V( x0)
!                if( V_trial < this%system%E_cut) then
!                    counter=counter+1
!                    Delta_E=0d0
!                    P = this%system%P(x0)
!                    U_move(k)= P
!                    !print*, "Umove prior to loop is ", U_move(k), "allocated?", allocated(U_move), " Valloc? ", allocated(this%Vij)
!            !  !$OMP parallel default(none) shared(Npoints,this%Vij,this%X,k,U_move,this%system%X,Natoms,x0) private(j,diff) reduction(+:Delta_E)
!            !  $OMP p#arallel default(none) shared(Npoints,Vij,X,k,U_move,x0) private(j,diff) reduction(+:Delta_E)
!            !$OMP parallel default(none) shared(Npoints,this,k,U_move,x0) private(j,diff) reduction(+:Delta_E)
!            !$OMP do 
!                    do j=1,Npoints
!                        !write(*,*) "k in parallel is", k, "j is", j
!                        if(j.ne.k) then
!                            U_move(j)=this%Pair_LJ_NRG( this%X(:,j), x0)
!                            Delta_E = Delta_E+ ( this%Vij(j,k) - U_move(j))
!                            !write(*,*) j,k,"Ujk=", this%Vij(j,k), " DE= ", Delta_E, " trial U_move was", U_move(j), &
!                            !    & "leading to ", this%Vij(j,k)-U_move(j)
!                        endif
!                    enddo
!            !$OMP end do
!            !$OMP end parallel
!                    !write(*,*) "Delta_E from reduce is", Delta_E
!                    !write(*,*) "new is", U_move(k), "old was", this%Vij(k,k), "P(x) is", P
!                    if(Delta_E>0d0 .or. exp(Delta_E) > .5) then
!                        this%Vij(:,k)=U_move(:)
!                        this%Vij(k,:)=U_move(:)
!                        accept=accept+1
!                        newmove = 1
!                        !write(*,*) k, Delta_E, "Accepted! Displacement was", &
!                        !    & sqrt(sum((mv_cutoff*(2*s-1))**2)), "V(x0)", V, this%X(:,k), x0! " DX ", sqrt(sum((this%X(:,k)-x0(:))**2))
!                        SDX = SDX + sqrt(sum((this%X(:,k)-x0(:))**2))
!                        deltae1=deltae1+Delta_E
!                        E = E - deltae1
!                        V = V + (V_trial - this%system%V( this%X(:,k) ))
!                        if (exp(Delta_E) > .5) then
!                            hotaccept=hotaccept+1
!                        endif
!                        this%X(:,k)=x0(:)
!                    endif
!                endif
!            enddo
        
        !if(mod(i,MMC_freq)==0)then
        !!character(len=*), parameter :: fmt = '(I10, "/", I10, " cLJ ", E6.3, " DX=", E6.3, " GE=", E7.4, " GDE=", E7.4," SE ",'&
        !!    &'E7.4, " <SDX> ", E6.3, " acc=", I7, "/", I8)'
        !    print fmt, i, N_MC, saved, this%cLJ, mv_cutoff, E, deltae1, V, SDX, accept-hotaccept, hotaccept, &
        !        & dble(accept + hotaccept)/dble(Npoints * MMC_freq)
        !    if(dble(accept)/counter<0.5)then 
        !        mv_cutoff=mv_cutoff*0.999
        !    else
        !        mv_cutoff=mv_cutoff*1.001
        !    endif
        !    mv_cutoff = max( mv_cutoff, 1e-2)
        !    !if (mv_cutoff < .0001) then
        !    !    mv_cutoff = .0001
        !    !endif
        !    !call Moments(Moment,x)
        !    !write(*,*) 'Moments:'
        !    !write(*,*) Moment(1:9)
        !    !write(18,*) plt_count, mv_cutoff
        !    !write(19,*) plt_count, deltae1
        !    if ( abs(deltae1) < this%opt_eps .and. dble(accept)/dble(Npoints*MMC_freq) > .2 ) then
        !        print*, "**CONVERGED ENERGY**"
        !        return
        !    end if
        !    !if ( mv_cutoff < 1e-5) then
        !    !    print*, "**CONVERGED TO STATIONARY POINT**"
        !    !    return
        !    !end if
        !    if ( deltae1 == 0.0 ) then
        !        this%cLJ = this%cLJ * .999
        !    else
        !        this%cLJ = this%cLJ * 1.001
        !    endif
        !    accept=0
        !    hotaccept=0
        !    counter=0
        !    plt_count=plt_count+1
        !    deltae1=0.0
        !endif
        saved='N'
        !write(*,*) "wrote xyz"
        !open(unit=21,file='grid.dat')
        !    do j=1,Npoints
        !        write(21,*) x(:,j)
        !    enddo
        !close(21)
    enddo
    !close(18)
    !close(19)

    end subroutine

    subroutine optimize_MC( this)
        class (grid_LJ_t) :: this
        type (grid_LJ_t)  :: obj

        double precision, dimension(:), allocatable :: s, x0
        double precision, dimension(:), allocatable :: U_move
        double precision :: Delta_E, mv_cutoff=1.0, deltae1=0.0, V, V_trial, P, E, SDX
        integer :: i, j, k, counter=0, Npoints, accept=0, hotaccept=0,plt_count=0, MMC_freq, N_MC, newmove=0, framesaved=1
        integer :: saveevery
        character(1) :: saved = 'N'
        !double precision, dimension(this%points_nr,this%points_nr)  :: Vij
        !double precision, dimension(3,this%points_nr)  :: X
        !character(len=*), parameter :: fmt = "(I8,A1,I10,A5,F6.3,A5,E8.6,A5,E8.6,A5,E8.6,5A,E8.6,A7,F6.3,5A,I4,I6)"
        character(len=*), parameter :: fmt = '(I10, " /", I10,1x, A1, "  cLJ=", ES13.6, "   DX=", ES13.6, "   GE=", ES13.6, &   
            &"   GDE=", ES13.6,"   SE=", ES13.6, "   <SDX>= ", ES13.6, " cold=", I4, "  hot=", I4, "   acc=", F8.4 )'

        !Vij = this%Vij
        !X   = this%X
        
        obj = this

        print*, "GRID_LJ OPTIMIZER", this%opt_tune, this%cLJ, " steps=", N_MC

        MMC_freq = this%opt_tune
        N_MC     = this%opt_steps
        !N_MC = 20000
        mv_cutoff= this%mv_cutoff
        saveevery = this%saveevery 

        !open(unit=18,file='mv_cut.dat')
        !open(unit=19,file='delE.dat')

        Npoints = this%points_nr

        allocate( U_move(this%points_nr), s(this%system%d), x0(this%system%d))

        open(unit=20,file='grid.xyz', access='APPEND')
        close(20)

        do i=1,Npoints
            do j=1, i-1
                this%Vij(j,i) = this%Pair_LJ_NRG( this%X(:,j), this%X(:,i))
            enddo
            do j=1, i-1
                this%Vij(i,j) = this%Vij(j,i)
            enddo
            this%Vij(i,i) = 0.0
        enddo
        i = 0
        do
            i = i + 1
            if ( N_MC == 0 .or. ( N_MC > -1 .and. i == N_MC)) then
                exit
            endif
            newmove = 0
            !write(*,*) i,N_MC
            E = sum(this%Vij) 
            V = 0.0
            SDX = 0.0
            do k=1, Npoints
                V = V + this%system%V( this%X(:,k))
            enddo
            do k=1, Npoints

                call random_number(s) 
                x0=this%X(:,k)+mv_cutoff*(2.0*s-1.0)
                !print*, "DELTA ", mv_cutoff*(2.0*s-1.0)
                !print*, "x0 is", x0
                V_trial = this%system%V( x0)
                if( V_trial < this%system%E_cut) then
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
                    !write(*,*) "Delta_E from reduce is", Delta_E
                    !write(*,*) "new is", U_move(k), "old was", this%Vij(k,k), "P(x) is", P
                    !if(Delta_E>0d0 .or. exp(Delta_E) > .99) then
                    if(Delta_E>0d0) then
                        this%Vij(:,k)=U_move(:)
                        this%Vij(k,:)=U_move(:)
                        accept=accept+1
                        newmove = 1
                        !write(*,*) k, Delta_E, "Accepted! Displacement was", &
                        !    & sqrt(sum((mv_cutoff*(2*s-1))**2)), "V(x0)", V, this%X(:,k), x0! " DX ", sqrt(sum((this%X(:,k)-x0(:))**2))
                        SDX = SDX + sqrt(sum((this%X(:,k)-x0(:))**2))
                        deltae1=deltae1+Delta_E
                        E = E - deltae1
                        V = V + (V_trial - this%system%V( this%X(:,k) ))
                        if (exp(Delta_E) > .99) then
                            hotaccept=hotaccept+1
                        endif
                        this%X(:,k)=x0(:)
                    endif
                endif
            enddo
        if( newmove > 0 .and. mod( i, saveevery) == 0 ) then
            open(unit=20,file='grid.xyz', access='APPEND')
                write(20,*) Npoints
                write(20,*) "step", i, "frame", framesaved
                do j=1,Npoints
                    write(20,*) "C     ", this%X(:,j)
                enddo
            close(20)
            newmove = 0
            framesaved = framesaved + 1
            saved='Y'
        end if
        
        if(mod(i,MMC_freq)==0)then
        !character(len=*), parameter :: fmt = '(I10, "/", I10, " cLJ ", E6.3, " DX=", E6.3, " GE=", E7.4, " GDE=", E7.4," SE ",'&
        !    &'E7.4, " <SDX> ", E6.3, " acc=", I7, "/", I8)'
            print fmt, i, N_MC, saved, this%cLJ, mv_cutoff, E, deltae1, V, SDX, accept-hotaccept, hotaccept, &
                & dble(accept + hotaccept)/dble(Npoints * MMC_freq)
            if(dble(accept)/counter<0.5)then 
                mv_cutoff=mv_cutoff*0.999
            else
                mv_cutoff=mv_cutoff*1.001
            endif
            mv_cutoff = max( mv_cutoff, 1e-2)
            !if (mv_cutoff < .0001) then
            !    mv_cutoff = .0001
            !endif
            !call Moments(Moment,x)
            !write(*,*) 'Moments:'
            !write(*,*) Moment(1:9)
            !write(18,*) plt_count, mv_cutoff
            !write(19,*) plt_count, deltae1
            if ( abs(deltae1) < this%opt_eps .and. dble(accept)/dble(Npoints*MMC_freq) > .2 ) then
                print*, "**CONVERGED ENERGY**"
                return
            end if
            !if ( mv_cutoff < 1e-5) then
            !    print*, "**CONVERGED TO STATIONARY POINT**"
            !    return
            !end if
            !if ( deltae1 == 0.0 ) then
            !    this%cLJ = this%cLJ * .999
            !else
            !    this%cLJ = this%cLJ * 1.001
            !endif
            accept=0
            hotaccept=0
            counter=0
            plt_count=plt_count+1
            deltae1=0.0
        endif
        saved='N'
        !write(*,*) "wrote xyz"
        !open(unit=21,file='grid.dat')
        !    do j=1,Npoints
        !        write(21,*) x(:,j)
        !    enddo
        !close(21)
    enddo
    !close(18)
    !close(19)

    end subroutine

    subroutine optimize_MCCG( this)
        class (grid_LJ_t) :: this

        double precision, dimension(:), allocatable :: s, x0
        double precision, dimension(:), allocatable :: U_move
        double precision, dimension(this%system%d,this%points_nr)  :: dUdx, totaldUdx
        double precision, dimension(this%points_nr)  :: buff
        double precision :: Delta_E, mv_cutoff=1.0, mv_cutoff_ref=1.0, deltae1=0.0, Eprev, Vprev, V, V_trial, P, E, SDX, dx=1.0
        double precision :: Vs, Vsprev, F, Fprev, time1, time2, krand
        integer :: i, j, k, counter=0, Npoints, accept=0, hotaccept=0,plt_count=0, MMC_freq, N_MC, newmove=0, framesaved=1
        integer :: saveevery
        logical :: converged = .false., adaptive = .true.
        character(1) :: saved = 'N'
        logical, parameter :: timing = .false.
        character(len=*), parameter :: fmt = '(I10, " /", I10,1x, I10, "  F=", ES15.8, "  DF=", ES15.8, "  E=", ES15.8, &
            "   DE=", ES15.8, "   DX=", ES15.8, "   V=", ES15.8, "   dV=", ES15.8,"   Vg=", ES15.8, "   dVg= ", ES15.8 )'

        print*, "GRID_LJ OPTIMIZER", this%opt_tune, this%cLJ, " steps=", N_MC

        MMC_freq = this%opt_tune
        N_MC     = this%opt_steps
        mv_cutoff= this%mv_cutoff
        mv_cutoff_ref= this%mv_cutoff
        saveevery = this%saveevery 

        Npoints = this%points_nr

        allocate( U_move(this%points_nr), s(this%system%d), x0(this%system%d))

        i = 0
        open(unit=20,file='grid.xyz', access='APPEND')
            write(20,*) Npoints
            write(20,*) "step", i, "frame", framesaved
            do j=1,Npoints
                write(20,*) "C     ", this%X(:,j)
            enddo
        close(20)
        !print*,"X allocated?", allocated( this%X)
        Eprev = 0.0
        Vprev = 0.0
        Vgprev = this%U()
        totaldUdx(:,:) = 0.0
        buff(:) = 0.0
        !$OMP parallel default(none) shared(this, Npoints) private(k)
        !$OMP do 
        do k=1,Npoints
            this%Vij(:,k) = 0.0
        enddo
        !$OMP end do
        !$OMP end parallel
        if( timing) print*, "SYSTEM V"
        call cpu_time( time1)
        !$OMP parallel default(none) shared(this, buff, Npoints) private(k)
        !$OMP do 
        do k=1,Npoints
            buff(k) = this%system%V( this%X(:,k))
            !do j=k+1, Npoints
            !    !print*, "X k j=", j, k, this%X(:,k), this%X(:,j)
            !    Vprev = Vprev + this%Pair_LJ_NRG( this%X(:,j), this%X(:,k))
            !    !print*, "V j k=", j, k, Vprev
            !enddo
        enddo
        !$OMP end do
        !$OMP end parallel

        do k=1,Npoints
            V = V + buff(k)
        enddo
        call cpu_time( time2)
        if( timing) print*, "TIME ", time2 - time1

        !$OMP parallel default(none) shared(this, buff, Npoints) private(k)
        !$OMP do 
        do k=1,Npoints
            buff(k) = sqrt(sum(this%dVdx( this%X(:,k))**2))
        enddo
        !$OMP end do
        !$OMP end parallel
        F = 0.0
        do j=1, Npoints
            F = F + buff(j)
        enddo
        F = F/Npoints

        Vg = this%U()
        dUdx(:,:) = 0.0
        !print*,"grid ", this%X
        do
            if ( N_MC == 0 .or. ( N_MC > -1 .and. i == N_MC)) then
                exit
            endif
            i = i + 1
            newmove = 0
            !print*, "eval grad", i
            if( timing) print*,"DVDX"
            call cpu_time( time1)

            ! THE MC PART
            call random_number(s)
            call random_number(krand)
            k = idnint(krand*Npoints)+1
            k = max(1, k)
            k = min(Npoints, k)
            x0=this%X(:,k)+mv_cutoff*(2.0*s-1.0)
            !print*, "DELTA ", mv_cutoff*(2.0*s-1.0)
            !print*, "x0 is", x0
            V_trial = this%system%V( x0)
            if( V_trial < this%system%E_cut) then
                counter=counter+1
                Delta_E=0d0
                P = this%system%P( x0)
                U_move(k) = P
                !print*, "Umove prior to loop is ", U_move(k), "allocated?", allocated(U_move), " Valloc? ", allocated(this%Vij)
        !  !$OMP parallel default(none) shared(Npoints,this%Vij,this%X,k,U_move,this%system%X,Natoms,x0) private(j,diff) reduction(+:Delta_E)
        !  $OMP p#arallel default(none) shared(Npoints,Vij,X,k,U_move,x0) private(j,diff) reduction(+:Delta_E)
                !$OMP parallel default(none) shared(Npoints,this,k,U_move,x0,V_trial) private(j,diff) reduction(+:Delta_E)
                !$OMP do 
                do j=1,Npoints
                    !write(*,*) "k in parallel is", k, "j is", j
                    if(j.ne.k) then
                        U_move(j)=this%Pair_LJ_NRG( this%X(:,j), x0)
                        Delta_E = Delta_E + ( U_move(j) - this%Vij(j,k) )
                        !write(*,*) j,k,"Ujk=", this%Vij(j,k), " DE= ", Delta_E, " trial U_move was", U_move(j), &
                        !    & "leading to ", this%Vij(j,k)-U_move(j)
                    endif
                enddo
                !$OMP end do
                !$OMP end parallel
                Delta_E = Delta_E !+ ( V_trial - this%system%V( this%X(:,k)) )
                !write(*,*) "Delta_E from reduce is", Delta_E
                !write(*,*) "new is", U_move(k), "old was", this%Vij(k,k), "P(x) is", P
                if(Delta_E<0d0 .or. exp(Delta_E) < 1.1) then
                !if( Delta_E<0d0) then
                    this%Vij(:,k)=U_move(:)
                    this%Vij(k,:)=U_move(:)
                    accept=accept+1
                    newmove = 1
                    Vg = Vg + Delta_E
                    V = V + ( V_trial - this%system%V( this%X(:,k)))
                    !F = F + ( sqrt( sum( this%dVdx( x0)**2)) - sqrt(sum( this%dVdx( this%X(:,k))**2)))
                    F = F + ( sqrt( sum( mv_cutoff*(2.0*s-1.0)**2)) - sqrt(sum( this%dVdx( this%X(:,k))**2)))
                    !write(*,*) k, Delta_E, "Accepted! Displacement was", &
                    !    & sqrt(sum((mv_cutoff*(2*s-1))**2)), "V(x0)", V, this%X(:,k), x0, Vg, Delta_E! " DX ", sqrt(sum((this%X(:,k)-x0(:))**2))
                    if (exp(Delta_E) < 1.1) then
                        hotaccept=hotaccept+1
                    endif
                    dUdx(:,k) = mv_cutoff*(2.0*s-1.0)
                else
                    !write(*,*) k, Delta_E, "REJECTED! Displacement was", &
                    !    & sqrt(sum((mv_cutoff*(2*s-1))**2)), "V(x0)", V, this%X(:,k), x0! " DX ", sqrt(sum((this%X(:,k)-x0(:))**2))
                    !dUdx(:,k) = 0.0
                endif
            endif

            this%X(:,k) = this%X(:,k) + dUdx(:,k)

            if( adaptive) then
            if( i > 0) then
                !if( sign(1.0d1,Eprev) .ne. sign(1.0d1,E) .or. V-Vprev > 0) then

                !if( V+Vg-Vprev-Vgprev > 0 ) then
                if( (dble(accept) / dble(i)) > .5 ) then
                    mv_cutoff = mv_cutoff * 1.1
                else
                    mv_cutoff = mv_cutoff * 0.9999
                endif
            endif
            else
                newmove=1
            endif

            E = mv_cutoff
            !print fmt,  i, N_MC, 'N', F, F - Fprev, V+Vg, V+Vg-Vprev-Vgprev, dble(accept)/dble(i), V, V-Vprev, Vg, Vg-Vgprev
            if( newmove > 0) then
                converged = .false.
                !if( abs(F-Fprev) < this%opt_eps) then
                !    print *, "force convergence" 
                !    converged = .true.
                !endif
                !if( abs(V+Vg - Vprev-Vgprev) < this%opt_eps) then
                !    print *, "energy convergence"
                !    converged = .true.
                !endif


                ! !$OMP parallel default(none) shared(this, buff, Npoints) private(k)
                ! !$OMP do 
                ! do k=1,Npoints
                !     buff(k) = sqrt(sum(this%dVdx( this%X(:,k))**2))
                !     !do j=k+1, Npoints
                !     !    !print*, "X k j=", j, k, this%X(:,k), this%X(:,j)
                !     !    Vprev = Vprev + this%Pair_LJ_NRG( this%X(:,j), this%X(:,k))
                !     !    !print*, "V j k=", j, k, Vprev
                !     !enddo
                ! enddo
                ! !$OMP end do
                ! !$OMP end parallel
                ! F = 0.0
                ! do j=1, Npoints
                !     F = F + buff(j)
                ! enddo
                ! F = F/Npoints



                if( mod( i, saveevery) == 0 .or. converged) then
                    print fmt,  i, N_MC, framesaved, F, F - Fprev, V+Vg, V+Vg-Vprev-Vgprev, E, V, V-Vprev, Vg, Vg-Vgprev
                    !print*, "step", i, " F= ", F, " dF= ", F - Fprev, " E= ", V+Vg, "DX", E, " V= ", V, " dV= ", V - Vprev, &
                    !    & " Vg= ", Vg, " dVg= ", Vg - Vgprev
                    open(unit=20,file='grid.xyz', access='APPEND')
                        write(20,*) Npoints
                        write(20,*) "step", i, "frame", framesaved, " E ", E, " V ", V, " Vg ", Vg, " F ", F
                        do j=1,Npoints
                            write(20,*) "C     ", this%X(:,j)
                        enddo
                    close(20)
                    saved='Y'
                    newmove = 1
                    framesaved = framesaved + 1
                endif
                if( converged) then
                    exit
                endif
                Eprev = E
                Vprev = V
                Vgprev = Vg
                Fprev = F
            else if (i == N_MC) then
                print fmt,  i, N_MC, -1, F, F - Fprev, V+Vg, V+Vg-Vprev-Vgprev, E, V, V-Vprev, Vg, Vg-Vgprev
            endif
            dUdx(:,k) = 0.0
        saved='N'
    enddo

    !call d2Vdx2(this, "hess.dat")
    !close(18)
    !close(19)
    !deallocate( U_move, s, x0)

    end subroutine

    subroutine d2Vdx2( this, fname)
        implicit none
        class (grid_LJ_t) :: this
        character(len=*) :: fname
        integer :: i,j,k,l
        integer :: N, d  
        real(8), dimension(this%system%d,2) :: ptA, ptB
        real(8), dimension(this%system%d) :: E
        real(8) :: x, y, a, b
        real(8), parameter :: dx = 1e-5
        real(8), dimension( this%points_nr*this%system%d, this%points_nr*this%system%d )  :: hessian
        integer, parameter :: fid = 16
        N = this%points_nr
        d = this%system%d

        print*,"Calculating hessian.."
        do i=1,N
            ptA(:,1) = this%X(:,i)
            ptA(:,2) = this%X(:,i)
            ptA(:,1) = this%dVdx( this%X(:,i))
            do j=1, i
                ptB(:,1) = this%X(:,j)
                ptB(:,2) = this%X(:,j)
                ptB(:,1) = this%dVdx( this%X(:,j))
                do k=1, d
                    x = ptA(k,2)
                    ptA(k,2) = ptA(k,2) + dx
                    do l=1, k
                        y = ptA(l,2)
                        ptA(l,2) = ptA(l,2) + dx
                        print*,"A", this%dVdx( ptA(:,2)), ptA(:,1) 
                        print*,"diff", (this%dVdx( ptA(:,2)) - ptA(:,1))
                        E = (this%dVdx( ptA(:,2)) - ptA(:,1)) / dx**2
                        hessian((i-1)*d + k, (j-1)*d + l) = E(k) * E(l)
                        hessian((j-1)*d + l, (i-1)*d + k) = hessian((i-1)*d + k, (j-1)*d + l)
                        ptA(l,2) = y
                    enddo
                    ptA(k,2) = x
                enddo
            enddo
        enddo

        print*,"Saving hessian to", fname
        open(fid, file=fname)
        write(fid,*) d, N
        do i=1, d*N
            write(fid,*) hessian(:,i)
        enddo
        close(fid)
    end subroutine

end module grid_LJ_class
