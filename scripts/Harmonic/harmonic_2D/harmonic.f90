!=============================================================================80
!                       Metropolis Monte Carlo Code Grid
!4-6-19 testinggggggggggggg
!I want to make a 1D and 2D version of this code.
!will be easier to start from the working 2D morse probably
!placing on github now and update in the future (4-26-19)
!take the 2D morse code and convert into 2d HO to better test eigenvalues
!this code seems to be fine (I generate a uniform grid and set alphas to be 
!a constant)
!==============================================================================!
!       Discussion:
!Quasi-Regular Gaussian Basis
!using optimized grid solve generalized eigenvalue problem
!==============================================================================!
!       Modified:
!   1 May 2019
!       Author:
!   Shane Flynn 
!==============================================================================!
module QRGB_mod
implicit none
!==============================================================================!
!                            Global Variables 
!==============================================================================!
!d              ==> Gaussian dsionality
!==============================================================================!
integer::d,NG
double precision,parameter::D_morse=12.
double precision,parameter::c_x=0.2041241
double precision,parameter::c_y=0.18371169
!==============================================================================!
contains
!==============================================================================!
function Potential(x)
!==============================================================================!
!       Discussion:
!Hard-coded Morse Potential Energy 
!==============================================================================!
implicit none
double precision::x(d),Potential
!Potential=D_morse*((exp(-c_x*x(1))-1)**2+(exp(-c_y*x(2))-1)**2)
Potential=0.5*x(1)**2+0.5*x(2)**2
end function Potential
!==============================================================================!
end module QRGB_mod
program main
use QRGB_mod
!==============================================================================!
!       Discussion:
!==============================================================================!
!d              ==> i-th gaussian dsionality (x^i=x^i_1,x^i_2,..,x^i_d)
!NG             ==> Number of basis functions
!x              ==>(d) all atom coordinates
!x0             ==>(d) store previous coordinate before trial move
!freq           ==> Interval to update mv_cutoff size
!accept         ==> number of accepted trial moves, for acceptance~50%  
!counter        ==> total number of moves, for acceptance~50%
!x2 ith gaussians coordinates for matrix elements
!z integration sequence form matrix elements
!x_ij= i-jth gaussian center
!z=sequence for integration
!t_i,t_f        ==> cpu time to ~ simulation time
!==============================================================================!
implicit none
character(len=50)::grid_in,alpha_in
integer::Nsobol,data_freq,n,i,j,counter,l
integer*8::skip
double precision::time1,time2,aij,r2
double precision,allocatable,dimension(:)::alpha,eigenvalues,x_ij
double precision,allocatable,dimension(:,:)::x,Smat,z,Vmat,Hmat
!==============================================================================!
!                       LLAPACK dsygv variables                                !
!==============================================================================!
integer::itype,info,lwork
double precision,allocatable,dimension(:)::work
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
call cpu_time(time1)
read(*,*) d
read(*,*) NG
read(*,*) grid_in
read(*,*) alpha_in
read(*,*) Nsobol
read(*,*) data_freq

Nsobol=Nsobol/data_freq
Nsobol=Nsobol*data_freq
write(*,*) 'data_freq=',data_freq, '  Nsobol=',Nsobol
skip=Nsobol
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(x(d,NG),x_ij(d),alpha(NG),eigenvalues(NG),Smat(NG,NG),Vmat(NG,NG))
allocate(z(d,Nsobol),Hmat(NG,NG))
write(*,*) 'Test 0; Successfully Read Input File'
!==============================================================================!
!                           Input GridPoints x(d,NG)
!==============================================================================!
open(16,File=grid_in)
do n=1,NG
    read(16,*) x(:,n)
enddo 
close(16)
!==============================================================================!
!                               Read Alphas
!==============================================================================!
open(17,File=alpha_in)
do n=1,NG
    read(17,*) alpha(n)
enddo 
close(17)
!==============================================================================!
!                           Overlap Matrix (S)
!==============================================================================!
do i=1,NG
    do j=i,NG
         aij=alpha(i)*alpha(j)/(alpha(i)+alpha(j))
         r2=sum((x(:,i)-x(:,j))**2)
         Smat(i,j)=(alpha(i)+alpha(j))**(-0.5*d)*exp(-0.5*aij*r2)
       Smat(j,i)=Smat(i,j)
    enddo
enddo
!==============================================================================!
!                   Check to see if S is positive definite
!If this is removed, you need to allocate llapack arrays before Hamiltonian 
!==============================================================================!
!lwork=-1
!lwork=3400
lwork=max(1,3*NG-1)
allocate(work(max(1,lwork)))
call dsyev('v','u',NG,Smat,NG,eigenvalues,work,Lwork,info)
!write(*,*) 'lwork', lwork
write(*,*) 'Info (Initial Overlap Matrix) ==> ', info
open(unit=19,file='overlap_eigenvalues.dat')
do i=1,NG
    write(19,*) eigenvalues(i)
enddo
close(19)
write(*,*) 'Test 4; Overlap Matrix is Positive Definite'
!==============================================================================!
!                   Generate Sequence For Evaluating Potential
!==============================================================================!
do l=1,Nsobol
    call sobol_stdnormal(d,skip,z(:,l))
enddo
write(*,*) 'Test 5; Successfully Generated Integration Sequence'
!==============================================================================! 
!                             Evaluate Potential
!==============================================================================!
Vmat=0d0
open(unit=21,file='eigenvalues.dat')
do counter=1,Nsobol/data_freq
   do i=1,NG
      do j=i,NG
         x_ij(:)=(alpha(i)*x(:,i)+alpha(j)*x(:,j))/(alpha(i)+alpha(j))
         do l=(counter-1)*data_freq+1,counter*data_freq
            Vmat(i,j)=Vmat(i,j)+Potential(x_ij(:)+z(:,l)/sqrt(alpha(i)+alpha(j))) 
         enddo
      enddo
   enddo
   !==============================================================================!
   !                   Solve Generalized Eigenvalue Problem
   !==============================================================================!
   write(*,*) 'Generalized Eigenvalue Iteration', counter*data_freq
   do i=1,NG
      do j=i,NG
         aij=alpha(i)*alpha(j)/(alpha(i)+alpha(j))
         r2=sum((x(:,i)-x(:,j))**2)
         Smat(i,j)=(alpha(i)+alpha(j))**(-0.5*d)*exp(-0.5*aij*r2)
         Smat(j,i)=Smat(i,j)
         ! kinetic energy:
         Hmat(i,j)=0.5*aij*(d-aij*r2)
         ! kinetic + potential energy
         Hmat(i,j)=(Hmat(i,j)+Vmat(i,j)/(counter*data_freq))*Smat(i,j)
         Hmat(j,i)=Hmat(i,j)
      enddo
   enddo
   itype=1
   eigenvalues=0d0
   call dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,eigenvalues,work,Lwork,info)
   write(21,*) eigenvalues(:)
   write(*,*) 'info ==> ', info
enddo
close(21)
!==============================================================================!
!                               output file                                    !
!==============================================================================!
call cpu_time(time2)
open(90,file='simulation.dat')
write(90,*) 'particle dsionality ==> ', d
write(90,*) 'NG ==> ', NG
write(90,*) 'Nsobol==>', Nsobol
write(90,*) 'Total Time ==> ', time2-time1
close(90)
write(*,*) 'Hello Universe!'
end program main
