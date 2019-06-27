!=============================================================================80
!                       Metropolis Monte Carlo Code Grid
!==============================================================================!
!       Discussion:
!Generate 2D uniform grid for testing purposes
!==============================================================================!
!       Modified:
!   30 April 2019
!       Author:
!   Shane Flynn 
!==============================================================================!
program main
implicit none
!==============================================================================!
integer::NG_1D,Nparticles,i,j,counter
double precision::lower,upper,alpha_par
double precision,allocatable,dimension(:)::points,alpha
double precision,allocatable,dimension(:,:)::x
lower=-10
upper=25
NG_1D=20
Nparticles=NG_1D**2
alpha_par=1d0
!==============================================================================!
allocate(points(NG_1D),x(2,Nparticles),alpha(Nparticles))
do i=1,NG_1D
    points(i)=lower+(i-1.)*(upper-lower)/(NG_1D-1.)
enddo
counter=1
do i=1,NG_1D
    do j=1,NG_1D
        x(1,counter)=points(i)
        x(2,counter)=points(j)
        counter=counter+1
    enddo
enddo
write(*,*) 'Test 3; Successfully Generated Gaussian Centers'
!==============================================================================!
!                       Generate Alpha Scaling 
!Input Paramater (same for every dimension)
!==============================================================================!
alpha=alpha_par
!==============================================================================!
!                   Write Gridpoints Configuration to file 
!==============================================================================!
open(unit=17,file='grid.dat')
do i=1,Nparticles
    write(17,*) x(:,i)
enddo
close(17)
write(*,*) 'Test 2; Successfully Generated Quasi-Regular Gridpoints'
!==============================================================================!
!                          Generate Alpha Scaling 
!==============================================================================!
open(unit=18,file='alphas.dat')
do i=1,Nparticles
    write(18,*) alpha(i)
enddo
close(18)
write(*,*) 'Test 3; Successfully Generated Gaussian Widths'
!==============================================================================!
write(*,*) 'Hello Universe!'
end program main
