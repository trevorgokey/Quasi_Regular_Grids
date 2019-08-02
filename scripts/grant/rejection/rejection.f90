!=============================================================================80
!                           Rejection  Method
!=============================================================================80
!               Discussion:
! Fortran 90 Implementation of the Rejection method for a pdf as described in
! fortran numerical methods. 
! A 2D Gaussian Potential, with uniform enclosing function. 
!               Modified:
! 17 September 2018
!               Author:
! Shane Flynn
!==============================================================================!
module rej_mod
implicit none
!==============================================================================!
!                               Global Variables 
!==============================================================================!
!               Integer:
! Dimen     ==> Dimensionality of the system
!==============================================================================!
integer ::  Dimen
contains
!==============================================================================!
!==============================================================================!
!==============================================================================!
subroutine pdf_eval(x,p)
!==============================================================================!
!       Discussion:
! 2D Gaussian Potential. [0,8], centered at 4, hard-coded (Adjust Accordingly)
!==============================================================================!
implicit none
double precision :: x(Dimen),p
p=0d0
p = exp(-(((x(1)-4)**2)/2+((x(2)-4)**2)/2))
!write(*,*) 'evaluate PDF', p
end subroutine pdf_eval
!==============================================================================!
end module rej_mod
!==============================================================================!
!==============================================================================!
!==============================================================================!
program reject
use rej_mod
implicit none
!==============================================================================!
!                               Variables 
!==============================================================================!
!               Integer:
! Nsobol        ==> Dimensionality of the system
! count_acc     ==> Number of Accepted Points
! count_rej     ==> Number of Rejections
!               Double Precision:
! pxy           ==> PDF Evaluation
! accept        ==> Randomly Generated Acceptance Condition 
! z             ==> (Dimen) Quasi-Random Points
! x             ==> (Dimen,Nsobol) Accepted Coordinates
!==============================================================================!
integer :: Nsobol,count_acc,count_rej,i
integer*8 :: ii,skip
double precision :: pxy,accept
double precision, allocatable :: z(:),x(:,:)
!==============================================================================!
!                           Read Input File
!==============================================================================!
read(*,*) Dimen
read(*,*) Nsobol
skip=Nsobol
allocate(z(Dimen),x(Dimen,Nsobol))
count_acc=0
count_rej=0
z=0d0
x=0d0
!==============================================================================!
!                   Store All Coordinates Generated
!==============================================================================!
open(unit=66,file='all_data.dat')
do while (count_acc < Nsobol)
    call sobol_stdnormal(Dimen,skip,z(:))
    write(66,*) z(:)
    call random_number(accept)
    call pdf_eval(z(:),pxy)
!==============================================================================!
!                       Accept/Reject Criteria 
!==============================================================================!
    if(pxy>accept) then
        count_acc = count_acc +1
        x(:,count_acc) = z(:)
    else
        count_rej = count_rej + 1
    endif
end do
!==============================================================================!
!                     Pipe Data to file (./a.out > data)
!==============================================================================!
do i=1,Nsobol
    write(*,*) x(:,i)
enddo
close(66)
!==============================================================================!
!                           Simulation Results 
!==============================================================================!
open(unit=67,file='output.dat')
write(67,*) 'Nsobol ', Nsobol
write(67,*) 'Dimensionality', dimen
write(67,*) 'Total Rejections', count_rej
close(67)
end program reject
