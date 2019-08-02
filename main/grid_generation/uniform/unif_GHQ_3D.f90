!=============================================================================80
!                      3D Morse Code with qLJ Basis 
!==============================================================================!
!       Discussion:
!use a truncated uniform grid to solve the 3D morse oscillator problem
!set alpha to be a constant use ~2000 grid points in total
! have the code make the uniform grid (very large) and then truncate points
!==============================================================================!
!       Modified:
!   15 May 2019
!       Author:
!   Shane Flynn 
!==============================================================================!
module QRGB_mod
implicit none
!==============================================================================!
!                            Global Variables 
!==============================================================================!
!d              ==> Gaussian Dimensionality
!==============================================================================!
double precision, parameter::pi = 4.*acos(-1d0)
integer::d,NG
double precision::E_cut,c_LJ,integral_P
!==============================================================================!
contains
!==============================================================================!
function Potential(x)
!==============================================================================!
!       Discussion:
!Hard-coded Morse Potential Energy 
!==============================================================================!
implicit none
double precision,parameter::D_morse=12.
double precision,parameter::c_x=0.2041241
double precision,parameter::c_y=0.18371169
double precision,parameter::c_z=0.16329928
double precision::x(d),Potential
Potential=D_morse*((exp(-c_x*x(1))-1)**2+(exp(-c_y*x(2))-1)**2+&
    (exp(-c_z*x(3))-1)**2)  
end function Potential
!==============================================================================!
function P_x(x)
!==============================================================================!
!       Discussion:
!Target Distribution Function
!P_x            ==> evaluate P(x)
!x_i            ==>(d) ith particles coordinate x^i_1,..,x^i_d
!Del_par        ==> Delta parameter for P(x) distribution, :=10% of Ecut
!integral_P     ==> Normalization factor for P(x)
!They set gamma=1 in the paper so I will just ignore it here
!==============================================================================!
implicit none 
double precision::Del_par,x(d),P_x
Del_par=0.01*E_cut
if(Potential(x)<E_cut) then
   P_x=(E_cut+Del_par-Potential(x))/integral_P
else        !set equal to 0 if beyond Ecut
   P_x=1d-8
end if
end function P_x
!==============================================================================!
subroutine normalize_P(N_Eval,xmin,xmax)
!==============================================================================!
!       Discussion:
!Normalize the target distribution function
!For Morse: integrating over [a,b],[a,b],[a,b] 
!int P(r)~Area_Square/N sum_n=1,N P(r_n)
!norm           ==> evaluate P(r)
!r              ==>(d) coordinates for evaluating P(r)
!N_eval         ==> number of evaluations for integral approximation
!a,b            ==> bounds for square to normalize over (b-a)^2 = area
!==============================================================================!
integer::N_eval,i
double precision::r(d),norm,xmin(d),xmax(d)
integral_P=1d0              !set equal to 1 so you can initially call P_x
norm=0d0                   !compute the normalization
do i=1,N_Eval      !only want points within the E_cut
   call random_number(r)
   r(:)=xmin(:)+r(:)*(xmax(:)-xmin(:))
   if(Potential(r)<E_cut) norm=norm+P_x(r)
enddo
do i=1,d
   norm=norm*(xmax(i)-xmin(i))
enddo
integral_P=norm/N_eval
end subroutine normalize_P 
!==============================================================================!
!==============================================================================!
end module QRGB_mod
!==============================================================================!
!==============================================================================!
program main
use QRGB_mod
!==============================================================================!
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
character(len=50)::grid_in,theory_in
integer::GH_order,n,i,j,l1,l2,l3,upper,lower,NGU,counter,k,NG_1D
double precision::time1,time2,aij,r2,alpha0,Vij
double precision,allocatable,dimension(:)::alpha,eigenvalues,x_ij,z,w,rr,points
double precision,allocatable,dimension(:,:)::x,Smat,Hmat,r,rsmall
double precision, allocatable :: xmin(:), xmax(:)
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
read(*,*) GH_order
read(*,*) alpha0
read(*,*) E_cut
read(*,*) c_LJ
read(*,*) grid_in
read(*,*) theory_in
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(x(d,NG),x_ij(d),rr(d),alpha(NG),eigenvalues(NG),Smat(NG,NG))
allocate(Hmat(NG,NG),z(GH_order),w(GH_order),xmin(d),xmax(d))
write(*,*) 'Test 0; Successfully Read Input File'
!==============================================================================!
!                           Compute Uniform Grid 
!==============================================================================!
NG_1D=50
NGU=NG_1D**3
upper=200
lower=-5
allocate(points(NG_1D),r(d,NGU),rsmall(d,NGU))
!==============================================================================!
!                       Generate Uniform Grid Points
!==============================================================================!
do i=1,NG_1D
    points(i)=lower+(i-1.)*(upper-lower)/(NG_1D-1.)
enddo
counter=1
do i=1,NG_1D
    do j=1,NG_1D
        do k=1,NG_1D
            r(1,counter)=points(i)
            r(2,counter)=points(j)
            r(3,counter)=points(k)
            counter=counter+1
        enddo
    enddo
enddo
open(unit=17,file='centers.dat')
do i=1,NGU
    write(17,*) r(:,i)
enddo
counter=0
rsmall=0d0
do i=1,NGU
    if(potential(r(:,i)).lt.E_cut)then
        counter=counter+1
        rsmall(:,counter)=r(:,counter)
    else
        write(*,*) 'potential',  potential(r(:,i)) 
    endif
enddo
write(*,*) 'NGU=', NGU
write(*,*) 'counter=', counter
open(unit=18,file='grid.dat')
do i=1,counter
    write(18,*) rsmall(:,i)
enddo
close(18)


























!!open(16,File=grid_in)
!!do n=1,NG
!!    read(16,*) x(:,n)
!!enddo 
!!close(16)
!!!==============================================================================!
!!!                   Determine Normalization constant for P_x
!!!==============================================================================!
!!xmin=minval(x,DIM=2)
!!xmax=maxval(x,DIM=2)
!!xmin(:)=xmin(:)-(xmax(:)-xmin(:))*0.1
!!xmax(:)=xmax(:)+(xmax(:)-xmin(:))*0.1
!!do i=1,d
!!   write(*,*) i,'  xmin,xmax=',xmin(i),xmax(i)
!!enddo
!!call normalize_P(10000,xmin,xmax)
!!==============================================================================!
!!                          Generate Alpha Scaling 
!!==============================================================================!
!!open(unit=17,file='alphas.dat')
!!do i=1,NG
!!    alpha(i)=alpha0/(c_LJ*(P_x(x(:,i))*NG)**(-1./d))**2
!!    write(17,*) alpha(i)
!!enddo
!!close(17)
!!write(*,*) 'Test 3; Successfully Generated Gaussian Widths'
!==============================================================================!
!                           Overlap Matrix (S)
!==============================================================================!
!!do i=1,NG
!!    do j=i,NG
!!         aij=alpha(i)*alpha(j)/(alpha(i)+alpha(j))
!!         r2=sum((x(:,i)-x(:,j))**2)
!!         Smat(i,j)=(sqrt(alpha(i)*alpha(j))/(alpha(i)+alpha(j)))**(0.5*d)&
!!             *exp(-0.5*aij*r2)
!!       Smat(j,i)=Smat(i,j)
!!    enddo
!!enddo
!!!==============================================================================!
!!!                   Check to see if S is positive definite
!!!If this is removed, you need to allocate llapack arrays before Hamiltonian 
!!!==============================================================================!
!!lwork=max(1,3*NG-1)
!!allocate(work(max(1,lwork)))
!!call dsyev('v','u',NG,Smat,NG,eigenvalues,work,Lwork,info)
!!write(*,*) 'Info (Initial Overlap Matrix) ==> ', info
!!write(*,*) 'RCN =', eigenvalues(1)/eigenvalues(NG)
!!open(unit=19,file='overlap_eigenvalues.dat')
!!do i=1,NG
!!    write(19,*) eigenvalues(i)
!!enddo
!!close(19)
!!!==============================================================================!
!!!             Use Gauss Hermit quadrature to evaluate the potential matrix
!!! z(GH-order) w(GH-order) --- quadrature points and weights
!!!==============================================================================!
!!call cgqf ( GH_order, 6, 0d0, 0d0, 0d0, 0.5d0, z, w )
!!w=w/sqrt(2*pi)
!!!==============================================================================!
!!!                   Solve Generalized Eigenvalue Problem
!!!==============================================================================!
!!do i=1,NG
!!  do j=i,NG
!!     aij=alpha(i)*alpha(j)/(alpha(i)+alpha(j))
!!     r2=sum((x(:,i)-x(:,j))**2)
!!     Smat(i,j)=(sqrt(alpha(i)*alpha(j))/(alpha(i)+alpha(j)))**(0.5*d)&
!!         *exp(-0.5*aij*r2)   
!!     Smat(j,i)=Smat(i,j)
!!     ! kinetic energy:
!!     Hmat(i,j)=0.5*aij*(d-aij*r2)
!!     ! potential energy Vij matrix element
!!     x_ij(:)=(alpha(i)*x(:,i)+alpha(j)*x(:,j))/(alpha(i)+alpha(j))
!!     Vij=0d0
!!     ! d=3
!!     do l1=1,GH_order
!!        do l2=1,GH_order
!!           do l3=1,GH_order
!!              rr(1)=z(l1)
!!              rr(2)=z(l2)
!!              rr(3)=z(l3)
!!              rr=x_ij+rr/sqrt(alpha(i)+alpha(j))
!!              Vij=Vij+w(l1)*w(l2)*w(l3)*Potential(rr)
!!           enddo
!!        enddo
!!     enddo
!!  ! kinetic + potential energy
!!     Hmat(i,j)=(Hmat(i,j)+Vij)*Smat(i,j)
!!     Hmat(j,i)=Hmat(i,j)
!!  enddo
!!enddo
!!itype=1
!!eigenvalues=0d0  ! VM ???
!!call dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,eigenvalues,work,Lwork,info)
!!write(*,*) 'info ==> ', info
!!open(unit=21,file='eigenvalues.dat')
!!write(21,*) alpha0, eigenvalues(:)
!!close(21)
!!open(unit=73,file='abs_error.dat')
!!open(unit=74,file='rel_error.dat')
!!open(25,File=theory_in)
!!do i=1,NG
!!   read(25,*) r2
!!   write(74,*) i, (eigenvalues(i)-r2)/r2
!!    write(73,*) i, abs(eigenvalues(i)-r2)
!!enddo
!!close(25)
!!close(73)
!!close(74)
!!!==============================================================================!
!!!                               output file                                    !
!!!==============================================================================!
!!call cpu_time(time2)
!!open(90,file='simulation.dat')
!!write(90,*) 'dimensionality ==> ', d
!!write(90,*) 'NG ==> ', NG
!!write(90,*) 'Integral_P==>', integral_P
!!write(90,*) 'c LJ==>', c_LJ
!!write(90,*) 'Ecut ==>', E_cut
!!write(90,*) 'alpha0==>', alpha0
!!do i=1,d
!!   write(90,*) i, ' xmin xmax ==>', xmin(i),xmax(i)
!!enddo
!!write(90,*) 'Total Time ==> ', time2-time1
!!close(90)
!!write(*,*) 'Hello Universe!'
end program main
