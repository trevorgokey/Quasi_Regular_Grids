!=============================================================================80
!                      2D Morse with Gauss Hermite Quadurature 
!Make 1 code for uniform grid and QLJ grid using quadrature
!we know the uniform grid code works, use this 1 code to test
!For now just read in uniform grid or read in qlh grid and alphas
!in the future when this is working need to compute alphas in this code itself
!==============================================================================!
!       Modified:
!   19 May 2019
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
double precision,parameter::pi=4.*atan(1d0)
double precision::E_cut,integral_P,c_LJ
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
double precision,parameter::omega(2)=(/0.2041241,0.18371169/)
double precision,parameter::D_morse=12.
Potential=D_morse*sum((exp(-omega(:)*x(:))-1.)**2 )
end function Potential
!==============================================================================!
function P(x)
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
double precision::x(d),P
double precision,parameter::Delta=0.1
if(Potential(x)<E_cut) then
    P=(E_cut*(1+Delta)-Potential(x))/integral_P          
else        !set equal to 0 if beyond Ecut
   P=1d-20
end if
end function P
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
   if(Potential(r)<E_cut) norm=norm+P(r)
enddo
do i=1,d
   norm=norm*(xmax(i)-xmin(i))
enddo
integral_P=norm/N_eval
end subroutine normalize_P 
!==============================================================================!
!==============================================================================!
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
character(len=50)::grid_in,theory_in
logical::unif_grid
integer:: GH_order,n,i,j,l1,l2
double precision::time1,time2,aij,r2,alpha0,Vij
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
call cpu_time(time1)
read(*,*) d
read(*,*) NG
read(*,*) GH_order   
read(*,*) alpha0
read(*,*) E_cut
read(*,*) c_LJ
read(*,*) grid_in
read(*,*) theory_in
read(*,*) unif_grid
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(x(d,NG),x_ij(d),rr(d),alpha(NG),eigenvalues(NG),Smat(NG,NG))
allocate(Hmat(NG,NG),theory(NG),z(GH_order),w(GH_order))
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
!                          Define Grid And Alphas
!==============================================================================!
if(unif_grid.EQV..TRUE.)then
    write(*,*) 'Using a Uniform Grid and Constant Alpha'
!==============================================================================!
!               Uniform Grid, no normalization, alpha=alpha0
!==============================================================================!
    alpha=alpha0
else
!==============================================================================!
!For now just read in alphas from sobol code and use the same normalization
!that worked in the sobol code
!==============================================================================!
    write(*,*) 'Generating Alphas for QLJ Grid'
    open(unit=17,file='alpha_sobol.dat')
    do i=1,NG
        read(17,*) alpha(i)
    enddo
    close(17)
    integral_P=2485.2238839034922  
!==============================================================================!
!In the future we will need to compute alphas and normalization in this code 
!==============================================================================!
!!!    xmin=minval(x,DIM=2)
!!!    xmax=maxval(x,DIM=2)
!!!    xmin(:)=xmin(:)-(xmax(:)-xmin(:))*0.1
!!!    xmax(:)=xmax(:)+(xmax(:)-xmin(:))*0.1
!!!    do i=1,d
!!!       write(*,*) i,'  xmin,xmax=',xmin(i),xmax(i)
!!!    enddo
!!!    call normalize_P(10000,xmin,xmax)
!!!    write(*,*) 'integral_P', Integral_P
!!!    do i=1,NG
!!!        alpha(i)=alpha0/(c_LJ*(P_x(x(:,i))*NG)**(-1./d))**2
!!!    enddo
endif
open(unit=17,file='alphas.dat')
!==============================================================================!
!                           Write Alphas to File
!==============================================================================!
open(unit=18,file='alphas_GHQ.dat')
do i=1,NG
    write(18,*) alpha(i)
enddo
close(18)
write(*,*) 'Test 3; Successfully Generated Alphas'
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
open(unit=19,file='overlap_eigenvalues.dat')
write(*,*) 'RCN =', eigenvalues(1)/eigenvalues(NG)
do i=1,NG
    write(19,*) eigenvalues(i)
enddo
close(19)
write(*,*) 'Test 4; Overlap Matrix is Positive Definite'
!==============================================================================!
!             Use Gauss Hermit quadrature to evaluate the potential matrix
! z(GH-order) w(GH-order) --- quadrature points and weights
!==============================================================================!
call cgqf ( GH_order, 6, 0d0, 0d0, 0d0, 0.5d0, z, w )
w=w/sqrt(2*pi)
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
     ! kinetic energy:
     Hmat(i,j)=0.5*aij*(d-aij*r2)
     ! potential energy Vij matrix element
     x_ij(:)=(alpha(i)*x(:,i)+alpha(j)*x(:,j))/(alpha(i)+alpha(j))
     Vij=0d0
     do l1=1,GH_order
        do l2=1,GH_order
           rr(1)=z(l1)
           rr(2)=z(l2)
           rr=x_ij+rr/sqrt(alpha(i)+alpha(j))
           Vij=Vij+w(l1)*w(l2)*Potential(rr) 
        enddo
     enddo
!==============================================================================!
!               In DGB formulaiton we normalize our Gaussians: P_ij(r)
!normalization factor is needed so scale V by (alpha(i)+alpha(j)/2pi)^{d/2}
!==============================================================================!
  ! kinetic + potential energy
     Hmat(i,j)=(Hmat(i,j)+Vij)*Smat(i,j)
     Hmat(j,i)=Hmat(i,j)
  enddo
enddo
itype=1
eigenvalues=0d0
call dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,eigenvalues,work,Lwork,info)
write(*,*) 'info ==> ', info
open(unit=21,file='eigenvalues.dat')
write(21,*) alpha0, eigenvalues(:)
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
write(73,*) alpha0, abs(theory(:)-eigenvalues(:))
close(73)
open(unit=74,file='rel_error.dat')
do i=1,NG
    write(74,*) i, (eigenvalues(i)-theory(i))/theory(i)
enddo
close(74)
!==============================================================================!
!                               output file                                    !
!==============================================================================!
call cpu_time(time2)
open(90,file='simulation.dat')
write(90,*) 'particle dsionality ==> ', d
write(90,*) 'NG ==> ', NG
write(90,*) 'alpha0==>', alpha0
write(90,*) 'Ecut==>', E_cut
write(90,*) 'GH_order==>', GH_order
write(90,*) 'Total Time ==> ', time2-time1
close(90)
write(*,*) 'Hello Universe!'
end program main
