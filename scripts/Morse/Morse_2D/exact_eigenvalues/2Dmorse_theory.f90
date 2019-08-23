!=============================================================================80
!                           2D Morse Eigenvalues
!=============================================================================80
!Compute Eigenvalues for the 2D Morse Oscillator (Analytic Expression)
!The eigenvalues will become negative due to the decreasing distance between 
!energy states at higher excitations (After 24 excitations in 1D). 
!==============================================================================!
!The energies output (theory.dat) contain all the combinatorics 
!These values are NOT sorted, see Python script "sort_2Dmorse.py" for sorting
!==============================================================================!
!       Modified:
!   5 May 2019
!       Author:
!   Shane Flynn
!==============================================================================!
program main
!==============================================================================!
!NB             ==>Number of Excitations (Number of Basis Functions)
!ax,ay          ==>Parameters Taken from JCP Garaschuk/Light 2001
!Dmorse         ==>Morse Parameter
!==============================================================================!
implicit none
integer::i,j,NB
double precision::D_morse,ax,ay
double precision,allocatable::eigx(:),eigy(:)
!==============================================================================!
!E(n):=a*sqrt(2*D_morse)*(n+1/2)-((a*sqrt(2*D_morse)*(n+1/2))^2/(4*D_morse))
!==============================================================================!
D_morse=12
ax=0.2041241  
ay=0.18371169 
NB=24                            !will become negative with more basis functions
!==============================================================================!
open(unit=17,file='1Dmorse_eig.dat')
allocate(eigx(NB),eigy(NB))
!==============================================================================!
!Compute 1D excitations, Excitaitons start at E_0 (ground state)
!==============================================================================!
do i=0,NB-1      
   eigx(i+1)=ax*sqrt(2.*D_morse)*(i+0.5)&
       - ((ax*sqrt(2.*D_morse)*(i+0.5))**2/(4.*D_morse))
   eigy(i+1)=ay*sqrt(2.*D_morse)*(i+0.5)&
       - ((ay*sqrt(2.*D_morse)*(i+0.5))**2/(4.*D_morse))
   write(17,*) eigx(i+1), eigy(i+1)
enddo
close(17)
!==============================================================================!
!Compute Combinatorics. Potential is Seperable i.e. E_total=E_x+E_y
!==============================================================================!
open(unit=90,file='2Dmorse_eig.dat')
do i=1,NB
    do j=1,NB
        write(90,*) eigx(i)+eigy(j)
    enddo
enddo
close(90)
write(*,*) 'Hello Universe!'
end program main
