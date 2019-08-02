!=============================================================================80
!                       Morse Eigenvalues
!=============================================================================80
!       Discussion:
!Compute Eigenvalues for the 3D Morse Oscillator (Exact Results)
!The eigenvalues will become negative due to the decreasing distance between
!energy states at higher excitations (After 24 excitations in 1D).
!E(n)=a*sqrt(2D)*(n+1/2)-((a*sqrt(2D)*(n+1/2))^2/(4D))
!==============================================================================!
!The energies output (theory.dat) contain all the combinatorics
!These values are NOT sorted, see "sort_3Dmorse.py" for sorting
!==============================================================================!
!       Modified:
!   15 May 2019
!       Author:
!   Shane Flynn
!==============================================================================!
program main
!==============================================================================!
!==============================================================================!
!==============================================================================!
!               Discussion:
!NB             ==> Number of Excitations
!omega          ==> Parameters Taken from JCP Garaschuk/Light.2001
!Dmorse         ==> Morse Parameter
!==============================================================================!
implicit none
integer::i,j,k,NB
double precision,parameter::D_morse=12d0
double precision,parameter::omega(3)=(/0.2041241,0.18371169,0.16329928/)
double precision,allocatable,dimension(:,:)::eigs
!==============================================================================!
NB=24
!==============================================================================!
open(unit=22,file='1Dmorse_eig.dat')
allocate(eigs(3,NB))
do i=0,NB-1
   eigs(:,i+1)=omega(:)*sqrt(2.*D_morse)*(i+0.5)-((omega(:)*sqrt(2.*D_morse)*&
       (i+0.5))**2/(4.*D_morse))
   write(22,*) eigs(:,i+1)
enddo
close(22)
!==============================================================================!
!Compute Combinatorics. Potential is Seperable i.e. E_total=E_x+E_y+E_z
!==============================================================================!
open(unit=23,file='3Dmorse_eig.dat')
do i=1,NB
    do j=1,NB
        do k=1,NB
            write(23,*) eigs(1,i)+eigs(2,j)+eigs(3,k)
        enddo
    enddo
enddo
close(23)
write(*,*) 'Hello Universe!'
end program main
