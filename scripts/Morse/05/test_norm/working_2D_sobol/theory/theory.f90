!=============================================================================80
!                       Morse Eigenvalues
!compute analytic eigenvalues for the 2D Morse Oscillator
!Shane Flynn
!working code 5-7-19 clean it up and place on github
!see https://ipfs.io/ipfs/QmXoypizjW3WknFiJnKLwHCnL72vedxjQkDDP1mXWo6uco/wiki/Morse_potential.html
!for the definition of the constants
!==============================================================================!
program main
!==============================================================================!
implicit none
integer::i,j,NB
double precision::D_morse,ax,ay
double precision,allocatable::eigx(:),eigy(:)
!==============================================================================!
!                       Theorectical Eigenvalues 2D Morse
!the first 24 are good, NB=24, after this 1D case becomes negative
!E(n)=a*sqrt(2D)*(n+1/2)-((a*sqrt(2D)*(n+1/2))^2/(4D))
!double precision,parameter::ax=0.2041241
!double precision,parameter::ay=0.18371169
!==============================================================================!
D_morse=12
ax=0.2041241  
ay=0.18371169 
NB=24
!==============================================================================!
open(unit=22,file='1D.dat')
allocate(eigx(NB),eigy(NB))
do i=0,NB
   eigx(i+1)=ax*sqrt(2.*D_morse)*(i+0.5) - ((ax*sqrt(2.*D_morse)*(i+0.5))**2/(4.*D_morse))
   eigy(i+1)=ay*sqrt(2.*D_morse)*(i+0.5) - ((ay*sqrt(2.*D_morse)*(i+0.5))**2/(4.*D_morse))
   write(22,*) eigx(i+1), eigy(i+1)
enddo
close(22)

open(unit=23,file='theory.dat')
do i=1,NB
    do j=1,NB
        write(23,*) eigx(i)+eigy(j)
    enddo
enddo
close(23)
write(*,*) 'Hello Universe!'
end program main
