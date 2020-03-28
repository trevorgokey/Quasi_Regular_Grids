
subroutine read_mrc14( fname)
    implicit none
    integer, parameter :: fid=15
    integer :: i
    character(len=*) :: fname
    integer :: NX, NY, NZ, MODE, NXSTART, NYSTART, NZSTART
    integer :: MX, MY, MZ, MAPC, MAPR, MAPS, MAP, MACHST, NLABL
    real(4) :: RMS, DMIN, DMAX, DMEAN
    real(4), dimension(3) :: CELLA, CELLB, CELLC, ORIGIN
    integer :: ISPG, NSYMBT 
    character(len=100):: EXTRA
    integer :: EXTRABUF(25)
    character, dimension(20,10) :: LABEL
    integer :: EXTTYP, NVERSION 
    character(len=1024) :: header
    character(len=*), parameter :: fmt = "A4"
    character,dimension(4) :: buf, buf2
    integer :: ibuf
    integer :: filesize
    real(4), dimension(:,:,:), allocatable :: voxel

    open(fid, file=fname, form='unformatted', access='direct', recl=4)

    read( unit=fid, rec=1)  NX
    read( unit=fid, rec=2)  NY
    read( unit=fid, rec=3)  NZ
    print*,"NX NY NZ",NX, NY, NZ
    read( unit=fid, rec=4)  MODE
    print*,"MODE", MODE
    read( unit=fid, rec=5)  NXSTART
    read( unit=fid, rec=6)  NYSTART
    read( unit=fid, rec=7)  NZSTART
    print*,"NXSTART NYSTART NZSTART",NXSTART, NYSTART, NZSTART
    read( unit=fid, rec=8)  MX
    read( unit=fid, rec=9)  MY
    read( unit=fid, rec=10) MZ
    print*,"MX MY MZ",MX, MY, MZ

    read( unit=fid, rec=11) ibuf
    CELLA(1) = transfer( ibuf, real(4))
    read( unit=fid, rec=12) ibuf
    CELLA(2) = transfer( ibuf, real(4))
    read( unit=fid, rec=13) ibuf
    CELLA(3) = transfer( ibuf, real(4))
    read( unit=fid, rec=14) ibuf
    CELLB(1) = transfer( ibuf, real(4))
    read( unit=fid, rec=15) ibuf
    CELLB(2) = transfer( ibuf, real(4))
    read( unit=fid, rec=16) ibuf
    CELLB(3) = transfer( ibuf, real(4))
    print*,"CELLA", CELLA
    print*,"CELLB", CELLB

    read( unit=fid, rec=17) MAPC
    read( unit=fid, rec=18) MAPR
    read( unit=fid, rec=19) MAPS
    print*,"MAP C/R/S", MAPC, MAPR, MAPS
    read( unit=fid, rec=20) DMIN
    read( unit=fid, rec=21) DMAX
    read( unit=fid, rec=22) DMEAN
    print*,"DMIN DMAX DMEAN", DMIN, DMAX, DMEAN
    read( unit=fid, rec=23) ISPG
    read( unit=fid, rec=24) NSYMBT
    print*,"ISPG NSYMBT", ISPG, NSYMBT
    do i=25,49
        read( unit=fid, rec=i) EXTRABUF(i-24)
    enddo
    !EXTRA = transfer( EXTRABUF, character, 100)
    EXTTYP = EXTRABUF(3)
    NVERSION = EXTRABUF(4)
    read( unit=fid, rec=50) ORIGIN(1)
    read( unit=fid, rec=51) ORIGIN(2)
    read( unit=fid, rec=52) ORIGIN(3)
    print*,"ORIGIN", ORIGIN

    !read( unit=fid, rec=53) MAP
    !MAP = transfer( ibuf, character, 4)

    read( unit=fid, rec=54) MACHST
    read( unit=fid, rec=55) RMS
    read( unit=fid, rec=55) NLABL

    call fseek( fid, 0, 2)
    filesize = ftell( fid)

    !print*,"reading",filesize - NX*NY*NZ*4, "bytes; filesize=", filesize
    if (MODE == 2) then
        allocate(voxel(NX,NY,NZ))
        do i=0,NX*NY*NZ-1
            read( unit=fid, rec=i+1+(filesize/4) - NX*NY*NZ) buf
            voxel( mod(i,NX)+1, mod(i/NX,NY)+1, mod(i/(NY*NZ), NZ)+1) = transfer( buf, real(4))
        enddo
    endif
    print*,"READ", NX*NY*NZ*4, "BYTES"
    close( fid)
end subroutine 

program test
    implicit none
    call read_mrc14( "test.map")
end program test
