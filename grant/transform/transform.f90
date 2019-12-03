!=============================================================================80
!                       Beasley Springer Moro
!=============================================================================80
!           Discussion:
! Fortran 90 implementation of the BSO algorithm for converting standard normal
! distribution from points distributed uniformly on the unit hypercube. 
! The code assumes points have already been generated and are read in from 
! a file (see sobol.f90 for generating uniformly distributed sobol points). 
!           Modified:
! 16 September 2018
!           Author:
! Shane Flynn
!=============================================================================80
module bsm_mod
implicit none
contains
function beasley_springer_moro(u) result(x)
implicit none
integer:: j
double precision :: r
double precision, dimension(:), intent(in) :: u
double precision, dimension(size(u)) :: x, y
double precision, parameter, dimension(0:3) :: a = (/ &
        2.50662823884, &
        -18.61500062529, &
        41.39119773534, &
        -25.44106049637 /)
double precision, parameter, DIMENSION(0:3) :: b = (/ &
        -8.47351093090, &
        23.08336743743, &
        -21.06224101826, &
        3.13082909833 /)
double precision, parameter, dimension(0:8) :: c = (/ &
        0.3374754822726147, &
        0.9761690190917186, &
        0.1607979714918209, &
        0.0276438810333863, &
        0.0038405729373609, &
        0.0003951896511919, &
        0.0000321767881768, &
        0.0000002888167364, &
        0.0000003960315187 /)
!==============================================================================!
y = u - 0.5D0
do j = 1, size(u)
    if (abs(y(j)) < 0.42) then
        r = y(j)*y(j)
        x(j) = y(j)*(((a(3)*r + a(2))*r + a(1))*r + a(0))/((((b(3)*r + b(2))*r + b(1))*r + b(0))*r + 1)
    else 
        if (y(j) > 0) then 
            r = log(-log(1-u(j)))
        elseif (y(j) < 0) then 
            r = LOG(-LOG(u(j)))
        endif 
        x(j) = c(0) + r*(c(1) + r*(c(2) + r*(c(3) + r*(c(4) + r*(c(5) + r*(c(6) + r*(c(7) + r*c(8))))))))
    if (y(j) < 0) then 
        x(j) = -x(j)
    endif
endif
enddo
end function beasley_springer_moro
end module bsm_mod
!==============================================================================!
!==============================================================================!
!==============================================================================!
program unif_norm
use bsm_mod
!==============================================================================!
!                                 Variables
!==============================================================================!
!           Integer:
! Npoints   ==> Number of Points 
! dimen     ==> Dimensionality (points are a vector)
!           Double Precision:
! a         ==> (Npoints,dimen), Uniform points, read from file
! x_std     ==> (Npoints,Dimen), Normally distributed points
!==============================================================================!
!==============================================================================!
implicit none
integer :: Npoints,dimen,i
character(len=50) :: in_data
double precision, allocatable :: a(:,:), x_std(:,:)
!==============================================================================!
!                           Read Input File
!==============================================================================!
read(*,*) Npoints
read(*,*) dimen 
read(*,*) in_data
allocate(a(Npoints, dimen),x_std(Npoints,dimen))
!==============================================================================!
!                               Read Uniform Points
!==============================================================================!
do i=1,Npoints
    random a(i,:)
enddo
!==============================================================================!
!                       Convert Uniform to Norm Distribution 
!==============================================================================!
x_std = 0d0
do i=1,Npoints
    x_std(i,:) = beasley_springer_moro(a(i,:))
    write(*,*) x_std(i,:)
end do
end program unif_norm
