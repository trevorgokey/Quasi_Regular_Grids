!=============================================================================80
!                   3D morse Grid Generation with Rejection
!==============================================================================!
!       Discussion:
!use rejection to generate 3D morse grid
!==============================================================================!
!       Modified:
!   5 June 2019
!       Author:
!   Shane Flynn 
!==============================================================================!
module rej_mod
implicit none
!==============================================================================!
!                            Global Variables 
!==============================================================================!
!d              ==> Particle Dimensionality
!E_cut          ==> Energy Cutoff for Particle Center
!integral_P     ==> Area under the curve for P
!x,y min/max     ==> Domain for P(x)
!==============================================================================!
integer::d,Npoints
double precision::E_cut,delta,integral_P 
!==============================================================================!
contains
!==============================================================================!
subroutine sobol_unif(d,skip,x_unif,xmin,xmax)
use sobol
implicit none
integer::i
double precision::s(d),xmin(d),xmax(d)
!> dimension
INTEGER(kind = 4), INTENT(IN) :: d
!> number of initial points to be skipped
INTEGER(kind = 8), INTENT(IN) :: skip   
!> return an array of doubles, standard normal
DOUBLE PRECISION, DIMENSION(d), INTENT(OUT) :: x_unif
x_unif=i8_sobol(int(d, 8), skip)
!scale uniform distribution to span domain
x_unif=xmin+x_unif*(xmax-xmin)
end subroutine sobol_unif
!==============================================================================!
function V(x)
!==============================================================================!
!       Discussion:
!Potential Energy (Hard-Coded 3D Morse)
!V            ==> evaluate V(x,y)
!x            ==>(d) ith particles coordinate x^i_1,..,x^i_d
!D_morse        ==> morse constant
!==============================================================================!
implicit none 
double precision::x(d),V
double precision,parameter::omega(3)=(/0.2041241,0.18371169,0.16329928/)
double precision,parameter::D_morse=12.
V=D_morse*sum( (exp(-omega(:)*x(:))-1.)**2 )
end function V
!==============================================================================!
function P(x)
!==============================================================================!
!Target Distribution Function
!==============================================================================!
implicit none 
double precision::x(d),P
if(V(x)<E_cut) then
   P=(E_cut-V(x)+delta)/integral_P          
else        !set equal to 0 if beyond Ecut
   P=1d-20
end if
end function P
!==============================================================================!
subroutine Moments_Reg(Moment,N,xmin,xmax)
!==============================================================================!
!The moments are computed using a regular square grid. 
!This is the most accurate and the best for 3d case
!==============================================================================!
integer :: N,i1,i2,i3,j
double precision :: r(d),xmin(d),xmax(d),Moment(0:9),dummy
Moment=0d0
do i1=0,N
   do i2=0,N
      do i3=0,N
         r(1)=xmin(1)+i1*(xmax(1)-xmin(1))/N
         r(2)=xmin(2)+i2*(xmax(2)-xmin(2))/N
         r(3)=xmin(3)+i3*(xmax(3)-xmin(3))/N
         dummy=P(r)
         Moment(0)=Moment(0)+dummy
         do j=1,d
            Moment(j)=Moment(j)+dummy*r(j)
            Moment(d+j)=Moment(d+j)+dummy*r(j)**2
         enddo
         Moment(7)=Moment(7)+dummy*r(1)*r(2)
         Moment(8)=Moment(8)+dummy*r(2)*r(3)
         Moment(9)=Moment(9)+dummy*r(3)*r(1)
      enddo
   enddo
enddo
dummy=1./N**3
do j=1,d
   dummy=dummy*(xmax(j)-xmin(j))
enddo
integral_P=dummy*Moment(0)
Moment(1:9)=Moment(1:9)/Moment(0)
end subroutine Moments_Reg
!==============================================================================!
subroutine Moments_MMC(Moment,N_MC,xmin,xmax)
!==============================================================================!
! Compute the moments using Metropolis MC and determin the box sizes (xmin,xmax)
!==============================================================================!
integer::N_MC,i,j
double precision::Moment(0:9),dummy,r_trial(d),r(d),s(d),xmin(d),xmax(d)
double precision,parameter::mv_cutoff=0.1
Moment=0d0
r=0d0
xmin=r
xmax=r
do i=1,N_MC
!==============================================================================!
!                   Generate coordinates for Random move trial
!random numbers generated (0,1), make it (-1,1) ==> s=2*s-1
!==============================================================================!
   call random_number(s) 
   r_trial=r+mv_cutoff*(2*s-1)
!==============================================================================!
!               Test to see if you should accept 
!==============================================================================!
   call random_number(dummy) 
   if(P(r_trial)/P(r).ge.dummy) then
      r=r_trial
      do j=1,d
         if(xmin(j)>r(j)) xmin(j)=r(j)
         if(xmax(j)<r(j)) xmax(j)=r(j)
      enddo
   endif
   do j=1,d
      Moment(j)=Moment(j)+r(j)
      Moment(d+j)=Moment(d+j)+r(j)**2
   enddo
   Moment(7)=Moment(7)+r(1)*r(2)
   Moment(8)=Moment(8)+r(2)*r(3)
   Moment(9)=Moment(9)+r(3)*r(1)
enddo
Moment=Moment/N_MC
end subroutine Moments_MMC
!==============================================================================!
subroutine Moments(Moment,x)
!==============================================================================!
!Compute the moments for a set of points
!==============================================================================!
integer::i,j
double precision::Moment(0:9),x(d,Npoints)
Moment=0d0
do i=1,Npoints
   do j=1,d
      Moment(j)=Moment(j)+x(j,i)
      Moment(d+j)=Moment(d+j)+x(j,i)**2
   enddo
   Moment(7)=Moment(7)+x(1,i)*x(2,i)
   Moment(8)=Moment(8)+x(2,i)*x(3,i)
   Moment(9)=Moment(9)+x(3,i)*x(1,i)
enddo
Moment(1:9)=Moment(1:9)/Npoints
end subroutine Moments
!==============================================================================!
end module rej_mod
!==============================================================================!
!==============================================================================!
program main
use rej_mod
!==============================================================================!
implicit none
integer::i,count_acc,count_rej,n_MC_moments,N_reg
integer*8::ii,skip
double precision::accept,pxy,moment(0:9),time1,time2
double precision,allocatable::z(:),x(:,:),xmin(:),xmax(:)
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
call cpu_time(time1)
read(*,*) d
read(*,*) Npoints
read(*,*) E_cut
read(*,*) delta
read(*,*) N_MC_moments  !number of MC moves using Monte carlo to compute the moments
read(*,*) N_reg     !number of grid points each dimension using regular grid to compute moments
skip=Npoints
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(z(d),x(d,Npoints),xmin(d),xmax(d))
write(*,*) 'Test 0; Successfully Read Input File'
!==============================================================================!
!                           Run some Monte Carlo
!==============================================================================!
integral_P=1d0              !set equal to 1 so you can initially call P
call Moments_MMC(Moment,N_MC_moments,xmin,xmax)
write(*,*) 'Moments from Metropolis MC using N_MC=',N_MC_moments
write(*,*) Moment(1:9)
do i=1,d
 write(*,*) i,' xmin xmax ==>', xmin(i),xmax(i)
enddo
xmin=xmin-(xmax-xmin)*0.01
xmax=xmax+(xmax-xmin)*0.01
call Moments_Reg(Moment,N_reg,xmin,xmax)
write(*,*) 'Moments using regular integration:'
write(*,*) Moment(1:9)
write(*,*) 'Integral of P =', Integral_P
!==============================================================================!
!Begin Rejection Methods
!==============================================================================!
count_acc=0
count_rej=0
z=0d0
x=0d0
open(unit=66,file='all_data.dat')
!==============================================================================!
!                     Generate Points Until total=Npoints
!==============================================================================!
do while (count_acc .lt. Npoints)
    call sobol_unif(d,skip,z(:),xmin,xmax)
    call random_number(accept)
!==============================================================================!
!                    Use PDF to Accept/Reject Coordinates
!==============================================================================!
    if(P(z)/(E_cut+delta).gt.accept) then   !use for sobol + rejection
!    if(P(z).gt.0.00000001) then            !use for uniform sobol with cutoff
        count_acc=count_acc+1
        x(:,count_acc)=z(:)
    else
        count_rej=count_rej+1
    endif
enddo
close(66)
open(unit=99,file='grid_sobol.dat')
do i=1,Npoints
    write(99,*) x(:,i)
enddo
close(99)
!==============================================================================!
!                               output file                                    !
!==============================================================================!
call cpu_time(time2)
open(90,file='simulation.dat')
write(90,*) 'particle dimensionality ==> ', d
write(90,*) 'Npoints ==> ', Npoints
write(90,*) 'P Normalization Iterations (N_MC) ==> ', N_MC_moments
write(90,*) 'P(x) Normalization ==> ', integral_P
write(90,*) 'E_cut ==> ', E_cut
write(90,*) 'N_MC_moments', N_MC_moments  
write(90,*) '# gridpoints regular moments', N_reg     
do i=1,3
 write(90,*) i,' xmin xmax ==>', xmin(i),xmax(i)
enddo
write(90,*) 'Total Time ==> ', time2-time1
write(90,*) 'Total Rejections', count_rej
close(90)
write(*,*) 'Hello Universe!'
end program main
