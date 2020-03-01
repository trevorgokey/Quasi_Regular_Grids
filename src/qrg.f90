
program main

    use grid_types, only: grid_t, grid_LJ_init
    use system_types, only: system_t, system_basis_gaussian_init

    integer i, fid, grid_potential, system_potential, points_nr, system_type, opt_type, system_input, grid_input
    character( len=1024) :: arg
    character( len=1024) :: input_fnm, grid_from
    class (grid_t), allocatable :: grid
    class (system_t), allocatable :: system
    namelist /ctrl/ grid_potential, system_potential, points_nr, system_type, opt_type, system_input, grid_input, grid_from

    i = 1
    input_fnm = ''
    do
        call get_command_argument( i, arg)
        if ( len_trim( arg) == 0 ) exit
        select case ( arg)
            case ('-i')
            i = i + 1
            call get_command_argument( i, arg)
            if( len_trim( arg) == 0 ) then
                print *, "-i needs a input filename"
                stop
            end if
            input_fnm = trim(arg)
        end select
        i = i + 1
    end do
    
    if ( len_trim( input_fnm) == 0 ) then
        print *, "Input file not specified with -i"
        stop
    end if

    fid = 11
    open( fid, file=input_fnm, iostat=i, action='read')
    if ( i == 0) then
        print *, "Opened file: ", trim(input_fnm)
    else
        print *, "Could not open file '", trim(input_fnm), "' error=", i
    end if

    read( unit=fid, nml=ctrl)
    print*, "Read grid_potential=",grid_potential, "system_type=",system_type, "opt_type=", opt_type

    select case ( grid_potential)
        case (1)
            grid = grid_LJ_init(fid, points_nr, grid_input, grid_from)
        case (2)
            !call grid_sobol( gptr, fid)
        case (3)
    end select
    select case ( system_potential)
        case (1)
            call system_basis_gaussian_init( system, fid)
        case (2)
            !call grid_sobol( gptr, fid)
        case (3)
    end select

    if (i.ne.0) then
        print*, "Could not init grid. Exiting"
        stop
    else
        print*, "Init grid success,   ret=", i
    end if

    call grid%grid_set_system( system)
    call grid%optimize()


    close( fid)

end program main
