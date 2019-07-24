!=============================================================================80
!               2D Henon-Heiles Gaussian Basis, Gauss Hermite Quadrature
!Code uses provided grid to construct a Gaussian Basis Set and solve the 
!Generalized Eigenvalue Problem (integrals are computed using GHQ) 
!This code defines alpha(i), the gaussian widths as a function of the original
!distribution function P(x).
!See GHQ_2D_HH.f90 for the nearest neighbor definition
!See HH_2D_grid.f90 for the code that generates the initial grid.
!==============================================================================!
!       Modified:
!   19 May 2019
!       Author:
!   Shane Flynn 
!==============================================================================!
module GHQ_mod
implicit none
!==============================================================================!
!                            Global Variables
!==============================================================================!
!       Discussion:
!d              ==>i-th gaussian dsionality (x^i=x^i_1,x^i_2,..,x^i_d)
!integral_P     ==>Normalization of P(x), computed in HH_2D_grid.f90
!E_cut          ==> Energy Cutoff Contour
!c_LJ           ==> parameter for LJ
!==============================================================================!
integer::d
double precision::integral_P,E_cut,c_LJ
!==============================================================================!
contains
!==============================================================================!
function V(x)
!==============================================================================!
!       Discussion:
!Potential Energy (Hard-Coded Henon-Heiles 2D)
!==============================================================================!
implicit none
double precision::x(d),V
double precision,parameter::lambda=sqrt(0.0125)
V=0.5*(x(1)**2+x(2)**2)+lambda*(x(1)**2*x(2)-x(2)**3/3.)
end function V
!==============================================================================!
function P(x)
!==============================================================================!
!       Discussion:
!Target Distribution Function
!==============================================================================!
!P              ==>evaluate P(x)
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!integral_P     ==>Normalization factor for distribution P(x)
!==============================================================================!
implicit none
double precision::x(d),P
if(V(x)<E_cut) then
   P=(E_cut-V(x))/integral_P
else        !set equal to 0 if beyond Ecut
   P=1d-20
end if
end function P
!==============================================================================!
end module GHQ_mod
!==============================================================================!
program main
use GHQ_mod
!==============================================================================!
!       Discussion:
!==============================================================================!
!NG             ==>Number of basis functions
!x              ==>(d) all atom coordinates
!x_ij           ==>i-jth gaussian center
!alpha0         ==>Scaling Parameter for the Gaussian Widths
!==============================================================================!
implicit none
character(len=50)::grid_in
integer::NG,GH_order,i,j,l1,l2
double precision::aij,r2,Vij,alpha0
double precision,parameter::pi=4.*atan(1d0)
double precision,allocatable,dimension(:)::alpha,eigenvalues,x_ij,theory,z,w,rr
double precision,allocatable,dimension(:,:)::x,Smat,Hmat
!==============================================================================!
!                       LLAPACK dsygv variables                                !
!==============================================================================!
integer::itype,info,lwork
double precision,allocatable,dimension(:)::work
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
read(*,*) d
read(*,*) NG
read(*,*) GH_order   
read(*,*) grid_in
read(*,*) integral_P
read(*,*) E_cut
read(*,*) c_LJ
read(*,*) alpha0
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(x(d,NG),x_ij(d),rr(d),alpha(NG),eigenvalues(NG),Smat(NG,NG))
allocate(Hmat(NG,NG),theory(NG),z(GH_order),w(GH_order))
write(*,*) 'Test 0; Successfully Read Input File'
!==============================================================================!
!                           Read GridPoints x(d,NG)
!==============================================================================!
open(17,File=grid_in)
do i=1,NG
    read(17,*) x(:,i)
enddo 
close(17)
!==============================================================================!
!                           Generate Alphas 
!generate Alpha(i) as a function of P(x_i)
!==============================================================================!
do i=1,NG
    alpha(i)=alpha0/(c_LJ*(P(x(:,i))*NG)**(-1./d))**2
enddo
write(*,*) 'Test 1; Successfully Generated Alphas for QLJ Grid'
!==============================================================================!
!                          Write Alphas to File 
!==============================================================================!
open(unit=18,file='alphas.dat')
do i=1,NG
    write(18,*) alpha(i)
enddo
close(18)
!==============================================================================!
!                           Overlap Matrix (S)
!==============================================================================!
do i=1,NG
    do j=i,NG
         aij=alpha(i)*alpha(j)/(alpha(i)+alpha(j))
         r2=sum((x(:,i)-x(:,j))**2)
         Smat(i,j)=(sqrt(alpha(i)*alpha(j))/(alpha(i)+alpha(j)))**(0.5*d)&
             *exp(-0.5*aij*r2)   
       Smat(j,i)=Smat(i,j)
    enddo
enddo
!==============================================================================!
!                   Check to see if S is positive definite
!define various arrays for using llapack, see llapack dsyev documentation
!==============================================================================!
lwork=max(1,3*NG-1)
allocate(work(max(1,lwork)))
call dsyev('v','u',NG,Smat,NG,eigenvalues,work,Lwork,info)
write(*,*) 'Info (Initial Overlap Matrix) ==> ', info
open(unit=19,file='overlap_eigenvalues.dat')
write(*,*) 'RCN =', eigenvalues(1)/eigenvalues(NG)
do i=1,NG
    write(19,*) eigenvalues(i)
enddo
close(19)
write(*,*) 'Test 2; Overlap Matrix is Positive Definite'
!==============================================================================!
!           Use Gauss Hermit quadrature to evaluate the potential matrix
!            z(GH-order) w(GH-order) --- quadrature points and weights
!                   See gen_hermite_rule.f90 for documentation
!==============================================================================!
call cgqf(GH_order,6,0d0,0d0,0d0,0.5d0,z,w)
w=w/sqrt(2.*pi)
!==============================================================================!
!                   Solve Generalized Eigenvalue Problem
!==============================================================================!
do i=1,NG
  do j=i,NG
     aij=alpha(i)*alpha(j)/(alpha(i)+alpha(j))
     r2=sum((x(:,i)-x(:,j))**2)
!==============================================================================!
!                                 Overlap
!==============================================================================!
     Smat(i,j)=(sqrt(alpha(i)*alpha(j))/(alpha(i)+alpha(j)))**(0.5*d)&
         *exp(-0.5*aij*r2)   
     Smat(j,i)=Smat(i,j)
!==============================================================================!
!                                 Kinetic
!==============================================================================!
     Hmat(i,j)=0.5*aij*(d-aij*r2)
!==============================================================================!
!                               Potential (Vij)
!==============================================================================!
     x_ij(:)=(alpha(i)*x(:,i)+alpha(j)*x(:,j))/(alpha(i)+alpha(j))
     Vij=0d0
     do l1=1,GH_order
        do l2=1,GH_order
           rr(1)=z(l1)
           rr(2)=z(l2)
           rr=x_ij+rr/sqrt(alpha(i)+alpha(j))
           Vij=Vij+w(l1)*w(l2)*V(rr) 
        enddo
     enddo
!==============================================================================!
!                   Hamiltonian=Kinetic+Potential (symmetric)
!==============================================================================!
     Hmat(i,j)=(Hmat(i,j)+Vij)*Smat(i,j)
     Hmat(j,i)=Hmat(i,j)
  enddo
enddo
!==============================================================================!
!                               Eigenvalues
!==============================================================================!
itype=1
call dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,eigenvalues,work,Lwork,info)
write(*,*) 'info ==> ', info
open(unit=20,file='eigenvalues.dat')
do i=1,NG
write(20,*) eigenvalues(i)
enddo
close(20)
!==============================================================================!
!                               output file                                    !
!==============================================================================!
open(90,file='simulation.dat')
write(90,*) 'particle dsionality ==> ', d
write(90,*) 'NG ==> ', NG
write(90,*) 'GH_order==>', GH_order
write(90,*) 'Integral_P ==>', Integral_P
write(90,*) 'Ecut ==>', E_cut
write(90,*) 'C_LJ ==>', c_LJ
write(90,*) 'alpha0 ==> ', alpha0
close(90)
write(*,*) 'Hello Universe!'
end program main
