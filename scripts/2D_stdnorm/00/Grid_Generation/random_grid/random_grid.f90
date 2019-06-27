!=============================================================================80
!                           Gaussian Random Grid                               !
!==============================================================================!
!       Discussion:
!Generate a d-dimensional Multivariate Gaussian Distribution of points using the
!Fortran random number generator and the Beasley-Springer-Moro transformation.
!==============================================================================!
!       Modified:                                                              
!   16 April 2019                                                              
!       Author:                                                                
!   Shane Flynn                                                                
!==============================================================================!
program main
!==============================================================================!
!               Discussion:
!==============================================================================!
!t_i,t_f        ==> cpu time 
!d              ==> dimensionality of the ith point (x^i=x^i_1,x^i_2,..,x^i_d)
!Npoints        ==> Total Number of Points
!x              ==>(d,Npoints) points coordinates
!==============================================================================!
implicit none
integer::d,Npoints,i
double precision::t_i,t_f
double precision,allocatable,dimension(:,:)::x
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
call cpu_time(t_i)
read(*,*) Npoints
read(*,*) d
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(x(d,Npoints))
!==============================================================================!
!                              Generate Points 
!==============================================================================!
open(unit=98,file='random.dat')
do i=1,Npoints
    call rand_stdnormal(d,x(:,i))
    write(98,*) x(:,i)
enddo
close(98)
!==============================================================================!
!                               Output File                                    !
!==============================================================================!
call cpu_time(t_f)
open(99,file='simulation.dat')
write(99,*) 'particle dimensionality ==> ', d
write(99,*) 'Npoints ==> ', Npoints
write(99,*) 'Total Time ==> ', t_f-t_i
close(99)
write(*,*) 'Hello Universe!'
end program main
