!=============================================================================80
!                    2D Morse with Gauss Hermite Quadrature 
!Code reads in grid and alphas and computes the associated eigenvalues
!Grid and alphas need to be generated elsewhere and provided here
!no Delta parameter
!needs teh gen_hermite_rule.f90 code for GHQ
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
integer::d
double precision::integral_P,E_cut,c_LJ
!==============================================================================!
contains
!==============================================================================!
function V(x)
!==============================================================================!
!       Discussion:
!Hard-coded Morse Potential Energy 
!==============================================================================!
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!V              ==>evaluate V(x)
!==============================================================================!
implicit none
double precision::x(d),V
double precision,parameter::omega(2)=(/0.2041241,0.18371169/)
double precision,parameter::D_morse=12.
V=D_morse*sum((exp(-omega(:)*x(:))-1.)**2 )
end function V
!==============================================================================!
function P(x)
!==============================================================================!
!       Discussion:
!Target Distribution Function
!==============================================================================!
!P              ==>evaluate P(x)
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!integral_P     ==>Normalization factor for P(x)
!Del_par        ==>Delta parameter for P(x) distribution, :=10% of Ecut
!==============================================================================!
implicit none
double precision::x(d),P
!double precision,parameter::Delta=0.1
if(V(x)<E_cut) then
!   P=(E_cut*(1+Delta)-V(x))/integral_P
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
!d              ==> i-th gaussian dsionality (x^i=x^i_1,x^i_2,..,x^i_d)
!NG             ==> Number of basis functions
!x              ==>(d) all atom coordinates
!x_ij           ==> i-jth gaussian center
!==============================================================================!
implicit none
character(len=50)::grid_in,theory_in
logical::unif_grid
integer::NG,GH_order,i,j,l1,l2
double precision::aij,r2,Vij,alpha0,RCN
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
read(*,*) theory_in
read(*,*) unif_grid
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
open(90,File=grid_in)
do i=1,NG
    read(90,*) x(:,i)
enddo 
close(90)
!==============================================================================!
!                           Generate Alphas alpha(NG)
!==============================================================================!
if(unif_grid.EQV..TRUE.)then
    alpha=alpha0
    write(*,*) 'Test 1; Uniform Grid, Alpha:=Constant'
else
    do i=1,NG
        alpha(i)=alpha0/(c_LJ*(P(x(:,i))*NG)**(-1./d))**2
    enddo
    write(*,*) 'Test 1; Successfully Generated Alphas for QLJ Grid'
endif
!==============================================================================!
!                          Write Alphas to File 
!==============================================================================!
open(unit=91,file='alphas.dat')
do i=1,NG
    write(91,*) alpha(i)
enddo
close(91)
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
!If this is removed, you need to allocate llapack arrays before Hamiltonian 
!==============================================================================!
lwork=max(1,3*NG-1)
allocate(work(max(1,lwork)))
call dsyev('v','u',NG,Smat,NG,eigenvalues,work,Lwork,info)
write(*,*) 'Info (Initial Overlap Matrix) ==> ', info
open(unit=92,file='overlap_eigenvalues.dat')
write(*,*) 'RCN =', eigenvalues(1)/eigenvalues(NG)
do i=1,NG
    write(92,*) eigenvalues(i)
enddo
close(92)
write(*,*) 'Test 2; Overlap Matrix is Positive Definite'
!==============================================================================!
!           Use Gauss Hermit quadrature to evaluate the potential matrix
!            z(GH-order) w(GH-order) --- quadrature points and weights
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
     Smat(i,j)=(sqrt(alpha(i)*alpha(j))/(alpha(i)+alpha(j)))**(0.5*d)&
         *exp(-0.5*aij*r2)   
     Smat(j,i)=Smat(i,j)
!kinetic energy:
     Hmat(i,j)=0.5*aij*(d-aij*r2)
!potential energy Vij matrix element
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
!kinetic + potential energy
     Hmat(i,j)=(Hmat(i,j)+Vij)*Smat(i,j)
     Hmat(j,i)=Hmat(i,j)
  enddo
enddo
itype=1
eigenvalues=0d0
call dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,eigenvalues,work,Lwork,info)
write(*,*) 'info ==> ', info
open(unit=21,file='eigenvalues.dat')
write(21,*) eigenvalues(:)
close(21)
!==============================================================================!
!                              Exact Eigenvalues
!==============================================================================!
open(25,File=theory_in)
do i=1,NG
    read(25,*) theory(i)
enddo 
close(25)
open(unit=73,file='abs_error.dat')
write(73,*) abs(theory(:)-eigenvalues(:))
close(73)
open(unit=74,file='rel_error.dat')
do i=1,NG
    write(74,*) i, (eigenvalues(i)-theory(i))/theory(i)
enddo
close(74)
!==============================================================================!
!                               output file                                    !
!==============================================================================!
open(90,file='simulation.dat')
write(90,*) 'particle dsionality ==> ', d
write(90,*) 'NG ==> ', NG
write(90,*) 'GH_order==>', GH_order
close(90)
write(*,*) 'Hello Universe!'
end program main
