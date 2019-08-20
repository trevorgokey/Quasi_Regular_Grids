!=============================================================================80
!                      Grid Generation via MMC Methods
!==============================================================================!
!       Discussion:
!Generate 2D Morse Grid and optimize with quasi-Lennard-Jones
!this code minimizes in an analagous manner to mmc methods
!see other versions for forced minimization
!keep this version
!rewrite poitential to use allocatabel vectors  
!need to properly edit this code, but it will be kept
!==============================================================================!
!       Modified:
!   30 April 2019
!       Author:
!   Shane Flynn 
!==============================================================================!
module qlj_mod
implicit none
!==============================================================================!
!                            Global Variables 
!==============================================================================!
!d              ==>Particle Dimensionality
!N_grid         ==>Number of gridpoints
!c_LJ           ==>parameter for LJ
!E_cut          ==>Energy Cutoff for Particle Center
!integral_P     ==>Area under the curve for P(x)
!x,y min/max    ==>Domain for normalizing P(x)
!==============================================================================!
integer::d,N_grid
double precision::c_LJ,E_cut,integral_P,xmin,xmax,ymin,ymax
!==============================================================================!
contains
!==============================================================================!
function random_integer(Nmin,Nmax) 
!==============================================================================!
!       Discussion:
!Randomly generate an integer in the range 1-N_grid
!       Variables:
!Nmin           ==>minimum index value (1)
!Nmax           ==>maximum index value (N_grid)
!random_integer ==>integer returned
!a              ==>Fortran intrinsic random number [0,1]
!==============================================================================!
implicit none
integer::Nmin,Nmax,random_integer
double precision::a
call random_number(a)
random_integer=floor(a*(Nmax-Nmin+1))+Nmin
end function random_integer
!==============================================================================!
function V(x)
!==============================================================================!
!       Discussion:
!Potential Energy (Hard-Coded 2D Morse)
!==============================================================================!
!V              ==>evaluate V(x,y)
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!D_morse        ==>morse constant
!==============================================================================!
implicit none 
double precision::x(d),V
double precision,parameter::omega_x=0.2041241
double precision,parameter::omega_y=0.18371169
double precision,parameter::D_morse=12.
V_x=D_morse*((exp(-omega_x*x_i(1))-1)**2+(exp(-omega_y*x_i(2))-1)**2)
end function V
!==============================================================================!
subroutine normalize_P(N_Eval)
!==============================================================================!
!       Discussion:
!Normalize the target distribution function
!integrating over the square [a,b],[a,b]: i.e. [-5,20], [-5,20]
!int P(r)~Area_Square/N sum_n=1,N P(r_n)
!==============================================================================!
!norm           ==> evaluate P(r)
!r              ==>(d) coordinates for evaluating P(r), sobol sequence
!N_eval         ==> number of evaluations for integral approximation
!a,b            ==> bounds for square to normalize over (b-a)^2 = area
!==============================================================================!
integer::N_eval,counter
double precision::r(2),norm
integral_P=1d0              !set equal to 1 so you can initially call P_x
norm=0d0                    !compute the normalization
counter=0
do while(counter.lt.N_Eval)      !only want points inside E_cut
    call random_number(r)
    r(1)=xmin+r(1)*(xmax-xmin)
    r(2)=xmin+r(2)*(ymax-ymin)
    if(V_x(r)<E_cut)then
        norm=norm+P_x(r)
        counter=counter+1
    endif
enddo
norm=norm*(xmax-xmin)*(ymax-ymin)/N_eval
integral_P=norm
end subroutine normalize_P 
!==============================================================================!
function P(x)
!==============================================================================!
!       Discussion:
!Target Distribution Function
!==============================================================================!
!P              ==> evaluate P(x)
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!Del_par        ==> Delta parameter for P(x) distribution, :=10% of Ecut
!integral_P     ==> Normalization factor for P(x)
!==============================================================================!
implicit none 
double precision::Del_par,x_i(d),P_x
Del_par=0.01*E_cut
if(V_x(x_i)<E_cut) then
   P_x=(E_cut+Del_par-V_x(x_i))/integral_P
else        !set equal to 0 if beyond Ecut
   P_x=1d-8
end if
end function P_x
!==============================================================================!
function Pair_LJ_NRG(x1,x2)
!==============================================================================!
!       Discussion:
!computes quasi LJ energy between 2 particles
!==============================================================================!
!x1             ==>(d) ith atoms coordinates
!x2             ==>(d) jth atoms coordinates
!a              ==>evaluate LJ
!sigma1/2       ==>c*sigma(P_x)
!Pair_LJ_NRG    ==>Energy of the i-j LJ potential
!==============================================================================!
implicit none 
double precision::x1(d),x2(d),a,b,Pair_LJ_NRG,sigma1,sigma2
a=sum((x1(:)-x2(:))**2)
sigma1=c_LJ*(P_x(x1)*N_grid)**(-1./d)    
sigma2=c_LJ*(P_x(x2)*N_grid)**(-1./d)    
b=(sigma2**2/a)**3
a=(sigma1**2/a)**3
Pair_LJ_NRG=a**2-a+b**2-b
end function Pair_LJ_NRG
!==============================================================================!
end module qlj_mod
!==============================================================================!
!==============================================================================!
program main
use qlj_mod
!==============================================================================!
!       Discussion:
!==============================================================================!
!N_Px           ==> Number of mmc iterations to normalize P_x
!NMC            ==> Number of MMC steps to execute
!x              ==>(dimen) all atom coordinates
!U              ==>(N_grid,N_grid) All i,j pair-wise energies (LJ-12-6)
!U_move         ==>(N_grid)  LJ-Energy associated with trial movement
!mv_cutoff      ==> maximum displacement parameter for trial displacement
!x0             ==>(d) store previous coordinate before trial move
!acc_coef       ==> scaling coefficient for making acceptance rate ~50%
!Delta_E        ==> change in energy due to potential
!t1             ==> random number for accepting higher energy movement
!s              ==>(d) random number to move coordinate by
!NMC_freq       ==> Interval to update mv_cutoff size
!accept         ==> number of accepted trial moves, for acceptance~50%  
!counter        ==> total number of moves, for acceptance~50%
!t_i,t_f        ==> cpu time to ~ simulation time
!==============================================================================!
implicit none
integer::NMC,NMC_freq,N_Px,accept,counter,i,j,k,n
double precision::Delta_E,mv_cutoff,t1,time1,time2,alpha0
double precision,allocatable,dimension(:)::x0,s,U_move,alpha
double precision,allocatable,dimension(:,:)::x,U
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
call cpu_time(time1)
read(*,*) d
read(*,*) N_grid
read(*,*) N_Px
read(*,*) NMC
read(*,*) NMC_freq
read(*,*) E_cut
read(*,*) c_LJ
read(*,*) alpha0
read(*,*) xmin
read(*,*) xmax
read(*,*) ymin
read(*,*) ymax
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(x(d,N_grid),x0(d),s(d),U(N_grid,N_grid),U_move(N_grid))
allocate(alpha(N_grid))
write(*,*) 'Test 0; Successfully Read Input File'
!==============================================================================!
!                   Determine Normalization constant for P_x
!==============================================================================!
call normalize_P(N_Px)
!==============================================================================!
!                       Generate Initial Distribution
!any point in the domain where Potential < Ecut
!==============================================================================!
n=0
do while(n.lt.N_grid)
    call random_number(s)
    s(1)=xmin+s(1)*(xmax-xmin)
    s(2)=xmin+s(2)*(ymax-ymin)
    if(V_x(s)<E_cut)then
        n=n+1
        x(:,n)=s(:)
    endif
enddo
open(unit=16,file='coor_ini.dat')
do i=1,N_grid
    write(16,*) x(:,i)
enddo
close(16)
write(*,*) 'Test 1; Successfully Generated Initial GridPoints'
!==============================================================================!
!                           Compute U[x_ij] 
!compute pairwise energies for all the initial GridPoints
!==============================================================================!
do i=2,N_grid
    do j=1,i-1
        U(i,j)=Pair_LJ_NRG(x(:,i),x(:,j))
        U(j,i)=U(i,j)
    enddo
enddo
!==============================================================================!
!                       Begin MMC to Optimize GridPoints
!==============================================================================!
accept=0
counter=0
mv_cutoff=0.01
do i=1,NMC
!==============================================================================!
!                       Randomly Select Atom to Move
!==============================================================================!
    k=random_integer(1,N_grid)
!==============================================================================!
!                   Generate coordinates for Random move trial
!random numbers generated (0,1), make it (-1,1) ==> s=2*s-1
!==============================================================================!
    call random_number(s) 
    x0=x(:,k)+mv_cutoff*(2*s-1)
!==============================================================================!
!                   Only consider point if V(trial) < Ecut 
!               Compute LJ Energy Change due to Trial Move
!==============================================================================!
    if(V(x0).lt.E_cut) then
        counter=counter+1
        U_move(k)=P_x(x0)
        Delta_E=0d0
        do j=1,N_grid
            if(j.ne.k) then
               U_move(j)=Pair_LJ_NRG(x(:,j),x0)
               Delta_E=Delta_E+U(j,k)-U_move(j)
            endif
        enddo
!==============================================================================!
!               Test to see if you should accept higher energy move
!if delta e > 0 than always larger than 1, always accept
!==============================================================================!
        call random_number(t1) 
        if(exp(Delta_E).ge.t1)then
           U(:,k)=U_move(:)
           U(k,:)=U_move(:)
           accept=accept+1
           x(:,k)=x0(:)
        endif
    endif
!==============================================================================!
!                           Update Cutoff Paramater
!acceptance rate ~50%, adjust random movement displacement length accordingly
!==============================================================================!
    if(mod(i,NMC_freq)==0)then
        write(*,*) 'MMC Iteration', n
        if(dble(accept)/counter<0.5)then 
            mv_cutoff=mv_cutoff*0.9
        else
            mv_cutoff=mv_cutoff*1.1
        endif
        accept=0
        counter=0
    endif
enddo
!==============================================================================!
!                   Write Gridpoints Configuration to file 
!==============================================================================!
open(unit=17,file='grid.dat')
do i=1,N_grid
    write(17,*) x(:,i)
enddo
close(17)
write(*,*) 'Test 2; Successfully Generated Quasi-Regular Gridpoints'
!==============================================================================!
!                          Generate Alpha Scaling 
!==============================================================================!
do i=1,N_grid
    alpha(i)=alpha0/(c_LJ*(P_x(x(:,i))*N_grid)**(-1./d))**2
enddo
open(unit=18,file='alphas.dat')
do i=1,N_grid
    write(18,*) alpha(i)
enddo
close(18)
write(*,*) 'Test 3; Successfully Generated Gaussian Widths'
!==============================================================================!
!                               output file                                    !
!==============================================================================!
call cpu_time(time2)
open(90,file='simulation.dat')
write(90,*) 'particle dimensionality ==> ', d
write(90,*) 'Number of gridpoints ==> ', N_grid
write(90,*) 'P_x Normalization Iterations (N_Px) ==> ', N_Px
write(90,*) 'P(x) Normalization ==> ', integral_P
write(90,*) 'Number of MMC Iterations for Grid ==> ', NMC
write(90,*) 'c_LJ ==> ', c_LJ
write(90,*) 'E_cut ==> ', E_cut
write(90,*) 'Final Move Cutoff ==> ', mv_cutoff
write(90,*) 'xmin ==>', xmin
write(90,*) 'xmax ==>', xmax
write(90,*) 'ymin ==>', ymin
write(90,*) 'ymax ==>', ymax
write(90,*) 'Total Time ==> ', time2-time1
close(90)
write(*,*) 'Hello Universe!'
end program main
