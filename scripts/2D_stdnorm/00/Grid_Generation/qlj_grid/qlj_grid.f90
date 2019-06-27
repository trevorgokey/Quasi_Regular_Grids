!=============================================================================80
!                       Metropolis Monte Carlo Grid
!==============================================================================!
!       Discussion:
!Generate a grid by treating grid points as particles and applying a modified 
!Lennard-Jones Potential for updating their position.
!Hard-Coded underlying distribution P(x), for demonstration purposes assume
!a gaussian distribution. 
!Initial distribution is generated using a gaussian-distributed sobol sequence
!Requires the Fortran sobol.f90 module
!==============================================================================!
!       Modified:
!   16 April 2019
!       Author:
!   Shane Flynn 
!==============================================================================!
module lj_mod
implicit none
!==============================================================================!
!                            Global Variables 
!==============================================================================!
!d              ==> Dimensionality of the ith point (x^i=x^i_1,x^i_2,..,x^i_d)
!epsilon_lj     ==> epsilon parameter for LJ 
!c_LJ           ==> parameter for LJ
!==============================================================================!
integer::d,Npoints
double precision::epsilon_LJ,c_LJ
double precision, parameter::pi=acos(-1d0)
!==============================================================================!
contains
!==============================================================================!
function random_integer(Nmin,Nmax) 
!==============================================================================!
!       Discussion:
!Randomly generate an integer in the range 1-Npoints
!       Variables:
!Nmin           ==> minimum index value (1)
!Nmax           ==> maximum index value (Npoints)
!random_integer ==> integer returned
!a              ==> Fortran intrinsic random number [0,1]
!==============================================================================!
implicit none
integer::Nmin,Nmax,random_integer
double precision::a
call random_number(a)
random_integer=floor(a*(Nmax-Nmin+1))+Nmin
end function random_integer
!==============================================================================!
function P_x(x)
!==============================================================================!
!       Discussion:
!computes the distribution for the i-th particle
!P_x            ==> evaluate P(x)
!x_i            ==>(d) ith particles coordinate x^i_1,..,x^i_d
!normalization factero
!==============================================================================!
implicit none 
double precision::x(d),P_x
!P_x=1.
P_x=(2.*pi)**(-d/2.)*exp(-0.5*sum(x(:)**2))
end function P_x
!==============================================================================!
function Pair_LJ_NRG(x1,x2)
!==============================================================================!
!       Discussion:
!computes LJ-12-6 energy between 2 particles
!       Variables:
!x_i            ==>(d) ith atoms coordinates
!x_j            ==>(d) jth atoms coordinates
!a              ==> evaluate LJ
!Pair_LJ_NRG    ==> Energy of the i-j LJ potential
!==============================================================================!
implicit none 
double precision::x1(d),x2(d),a,b,Pair_LJ_NRG,sigma1,sigma2
a=sum((x1(:)-x2(:))**2)
sigma1=c_LJ*(P_x(x1)*Npoints)**(-1./d)    
sigma2=c_LJ*(P_x(x2)*Npoints)**(-1./d)    
b=(sigma2**2/a)**3
a=(sigma1**2/a)**3
Pair_LJ_NRG=4.*epsilon_LJ*(a**2-a+b**2-b)
end function Pair_LJ_NRG
!==============================================================================!
end module lj_mod
!==============================================================================!
!==============================================================================!
program main
use lj_mod
!==============================================================================!
!               Discussion:
!==============================================================================!
!t_i,t_f        ==> cpu time to ~ simulation time
!Npoints     ==> Number of Particles
!N_iter         ==> Number of MMC steps to execute
!beta           ==> 1/kT, inverse temperature
!x              ==>(dimen) all atom coordinates
!dimen          ==> total system dimensionality (configuration space)
!d              ==> i-th particle dimensionality (x^i=x^i_1,x^i_2,..,x^i_d)
!V              ==>(Npoints) All energies due to distribution V(x)=-ln[P(x)]
!U              ==>(Npoints,Npoints) All i,j pair-wise energies (LJ-12-6)
!U_move         ==>(Npoints)  LJ-Energy associated with trial movement
!V_move         ==> V(x) associated with trial move V(x) = -ln[P(x)]
!mv_cutoff      ==> maximum displacement parameter for trial displacement
!x0             ==>(d) store previous coordinate before trial move
!acc_coef       ==> scaling coefficient for making acceptance rate ~50%
!max_move       ==> cutoff for local movement
!Delta_E        ==> change in energy due to potential
!t1             ==> random number for accepting higher energy movement
!s              ==>(d) random number to move coordinate by
!freq           ==> Interval to update mv_cutoff size
!accept         ==> number of accepted trial moves, for acceptance~50%  
!counter        ==> total number of moves, for acceptance~50%
!moment_i/f     ==> the second moment of the initial and final distributions
!ex_mom the exact moment determine from sobol sequence simulation
!==============================================================================!
implicit none
integer::i,j,k,n,N_iter,accept,counter,freq
integer*8::ii,skip                          !need *8 for sobol generator!
double precision::t1,Delta_E,mv_cutoff,beta,t_i,t_f,moment_i,moment_f,ex_mom,error
double precision,allocatable,dimension(:) :: x0,s,U_move
double precision,allocatable,dimension(:,:) :: x,U
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
call cpu_time(t_i)
read(*,*) Npoints
read(*,*) d
read(*,*) N_iter
read(*,*) beta
read(*,*) epsilon_LJ
read(*,*) c_LJ
read(*,*) mv_cutoff
read(*,*) freq
skip=Npoints
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(x(d,Npoints),x0(d),s(d),U(Npoints,Npoints),U_move(Npoints))
!==============================================================================!
!                           Generate Initial State
!use multivariant gaussian distribution (quasi random sequence)  to generate
!initial sequence, assume Nparticels*d configuration space (dimen)
!i.e.   x=([x_1,x_2,..,x_d],[x_1,x_2,..,x_d],..[x_1,x_2,..,x_d])
!==============================================================================!
open(unit=16,file='coor_ini.dat')
do ii=1,Npoints
    call sobol_stdnormal(d,skip,x(:,ii))
    write(16,*) x(:,ii)
enddo
close(16)
!==============================================================================!
!               Compute second moment of sobol distribution   <|x|2> 
!==============================================================================!
moment_i=0d0
do i=1,Npoints
    moment_i=moment_i+sum(x(:,i)**2)
enddo
moment_i=moment_i/Npoints
write(*,*) 'moment ', moment_i
!==============================================================================!
!                               Compute U[x_ij] 
!compute pairwise energies for the entire set of initial particle 
!==============================================================================!
do i=2,Npoints
    do j=1,i-1
        U(i,j)=Pair_LJ_NRG(x(:,i),x(:,j))
        U(j,i)=U(i,j)
    enddo
enddo
!==============================================================================!
!                               Compute V[x_i] 
!compute pairwise energies for the entire set of initial particle 
!==============================================================================!
do i=1,Npoints
    U(i,i)=-log(P_x(x(:,i)))
enddo
!==============================================================================!
!                               Begin MMC 
!==============================================================================!
accept=0
counter=0
do n=1,N_iter
!==============================================================================!
!                       Randomly Select Atom to Move
!==============================================================================!
    k=random_integer(1,Npoints)
!==============================================================================!
!                   Generate coordinates for Random move trial
!random numbers generated [0,1], make it [-1,1] ==> s=2*s-1
!==============================================================================!
    call random_number(s) 
    counter=counter+1
!==============================================================================!
!                           make trial move
!==============================================================================!
    x0=x(:,k)+mv_cutoff*(2*s-1)
!==============================================================================!
!                   Compute Energy Change due to Trial Move
!first compute Delta V = Vnew[x_i] - V[x_i]
!Then compute delta U and add together for total change in energy
!==============================================================================!
    U_move(k)=-log(P_x(x0))
    Delta_E=U(K,k)-U_move(k) 
    do j=1,Npoints
        if(j.ne.k) then
           U_move(j)=Pair_LJ_NRG(x(:,j),x0)
           Delta_E=Delta_E+U(j,k)-U_move(j)
        endif
    enddo
!==============================================================================!
!               Test to see if you should accept higher energy move
!criteria depends on how you define inequality
!if beta > 0 always, if delta e > 0 than always larger than 1, always accept
!==============================================================================!
    call random_number(t1) 
    if(exp(beta*Delta_E).ge.t1)then
       U(:,k)=U_move(:)
       U(k,:)=U_move(:)
       accept=accept+1
       x(:,k)=x0(:)
    else
!==============================================================================!
!                   Otherwise Reject Configuration
!==============================================================================!
    endif
!==============================================================================!
!                           Update Cutoff Paramater
!acceptance rate ~50%, adjust random movement displacement length accordingly
!==============================================================================!
    if(mod(n,freq)==0)then
        write(*,*) 'iteration', n
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
!                   Write Final Configuration to file 
!may want to make this as a function of MMC Iteration in the future
!==============================================================================!
open(unit=17,file='qlj_grid.dat')
do i=1,Npoints
    write(17,*) x(:,i)
enddo
!==============================================================================!
!                               output file                                    !
!==============================================================================!
call cpu_time(t_f)
open(90,file='simulation.dat')
write(90,*) 'Npoints ==> ', Npoints
write(90,*) 'particle dimensionality ==> ', d
write(90,*) 'Number of MMC Iterations ==> ', N_iter
write(90,*) 'epsilon LJ ==> ', epsilon_LJ 
write(90,*) 'c_LJ ==> ', c_LJ
write(90,*) 'Final Move Cutoff ==> ', mv_cutoff
write(90,*) 'Moment Sobol Dist.', moment_i
write(90,*) 'Moment final Dist.', moment_f
write(90,*) 'Moment Error', error
write(90,*) 'Total Time ==> ', t_f-t_i
close(90)
write(*,*) 'Hello Universe!'
end program main
