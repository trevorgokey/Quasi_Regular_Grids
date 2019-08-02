!=============================================================================80
!                       Metropolis Monte Carlo Code Grid
!==============================================================================!
!       Discussion:
!Compute distribion and eigenvalues for 2D morse as described in Garaschuk
!       TODO:
!==============================================================================!
!normalize PX                                   !DONE
!Generate Initial Distribution                  !DONE
!Remove single particle potential,epsilon,beta  !DONE
!==============================================================================!
!       Modified:
!   30 April 2019
!       Author:
!   Shane Flynn 
!==============================================================================!
module lj_mod
implicit none
!==============================================================================!
!                            Global Variables 
!==============================================================================!
!d              ==> Particle Dimensionality
!epsilon_lj     ==> epsilon parameter for LJ 
!c_LJ           ==> parameter for LJ
!E_cut          ==> Energy Cutoff for Particle Center
!D_morse        ==> Morse parameter
!integral_P     ==> Area under the curve for P_x
!x,ymin/max     ==> Domain for P(x)
!==============================================================================!
integer::d,Nparticles
double precision::c_LJ,E_cut,Del_par,integral_P,xmin,xmax,ymin,ymax
!==============================================================================!
contains
!==============================================================================!
function random_integer(Nmin,Nmax) 
!==============================================================================!
!       Discussion:
!Randomly generate an integer in the range 1-Nparticles
!       Variables:
!Nmin           ==> minimum index value (1)
!Nmax           ==> maximum index value (Nparticles)
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
function V_x(x_i)
!==============================================================================!
!       Discussion:
!Potential Energy (Hard-Coded 2D Morse)
!V_x            ==> evaluate V(x,y)
!x_i            ==>(d) ith particles coordinate x^i_1,..,x^i_d
!D_morse        ==> morse constant
!==============================================================================!
implicit none 
double precision::x_i(d),V_x
double precision,parameter::D_morse=12.
double precision,parameter::omega_x=0.2041241
double precision,parameter::omega_y=0.18371169
V_x=D_morse*((exp(-omega_x*x_i(1))-1)**2+(exp(-omega_y*x_i(2))-1)**2)
end function V_x
!==============================================================================!
subroutine normalize_P(N_Eval)
!==============================================================================!
!       Discussion:
!Normalize the target distribution function
!integrating over the square [a,b],[a,b]: i.e. [-5,20], [-5,20]
!int P(r)~Area_Square/N sum_n=1,N P(r_n)
!norm_P         ==> evaluate P(r)
!r              ==>(d) coordinates for evaluating P(r), sobol sequence
!N_eval         ==> number of evaluations for integral approximation
!a,b            ==> bounds for square to normalize over (b-a)^2 = area
!==============================================================================!
integer::N_eval,counter
double precision::r(2),norm
integral_P=1d0              !set equal to 1 so you can initially call P_x
norm=0d0                    !compute the normalization
counter=0
do while(counter.lt.N_Eval)      !only want points within the E_cut
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
function P_x(x_i)
!==============================================================================!
!       Discussion:
!Target Distribution Function
!P_x            ==> evaluate P(x)
!x_i            ==>(d) ith particles coordinate x^i_1,..,x^i_d
!Del_par        ==> Delta parameter for P(x) distribution
!integral_P     ==> Normalization factor for P(x)
!They set gamma=1 in the paper so I will just ignore it here
!==============================================================================!
implicit none 
double precision::x_i(d),P_x
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
!       Variables:
!x_i            ==>(d) ith atoms coordinates
!x_j            ==>(d) jth atoms coordinates
!a              ==> evaluate LJ
!Pair_LJ_NRG    ==> Energy of the i-j LJ potential
!==============================================================================!
implicit none 
double precision::x1(d),x2(d),a,b,Pair_LJ_NRG,sigma1,sigma2
a=sum((x1(:)-x2(:))**2)
sigma1=c_LJ*(P_x(x1)*Nparticles)**(-1./d)    
sigma2=c_LJ*(P_x(x2)*Nparticles)**(-1./d)    
b=(sigma2**2/a)**3
a=(sigma1**2/a)**3
Pair_LJ_NRG=a**2-a+b**2-b
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
!Nparticles     ==> Number of Particles
!N_iter         ==> Number of MMC steps to execute
!beta           ==> 1/kT, inverse temperature
!x              ==>(dimen) all atom coordinates
!dimen          ==> total system dimensionality (configuration space)
!d              ==> i-th particle dimensionality (x^i=x^i_1,x^i_2,..,x^i_d)
!V              ==>(Nparticles) All energies due to distribution V(x)=-ln[P(x)]
!U              ==>(Nparticles,Nparticles) All i,j pair-wise energies (LJ-12-6)
!U_move         ==>(Nparticles)  LJ-Energy associated with trial movement
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
!moment_i/f     ==> the moments of the distributions
!==============================================================================!
implicit none
integer::i,j,k,n,N_iter,accept,counter,freq,Nsobol,NMC
double precision::t1,Delta_E,mv_cutoff,t_i,t_f
double precision,allocatable,dimension(:)::x0,s,U_move
double precision,allocatable,dimension(:,:)::x,U
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
call cpu_time(t_i)
read(*,*) Nparticles
read(*,*) d
read(*,*) NMC
read(*,*) N_iter
read(*,*) c_LJ
read(*,*) E_cut
read(*,*) Del_par
read(*,*) mv_cutoff
read(*,*) freq
read(*,*) xmin
read(*,*) xmax
read(*,*) ymin
read(*,*) ymax
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(x(d,Nparticles),x0(d),s(d),U(Nparticles,Nparticles),U_move(Nparticles))
!==============================================================================!
!                   Determine Normalization constant for P_x
!==============================================================================!
write(*,*) 'Integral_P', integral_P
call normalize_P(NMC)
write(*,*) 'Integral_P', integral_P
!==============================================================================!
!                   Generate Initial Distribution
!generate a point in domain, if Potential is within Ecut keep 
!==============================================================================!
n=0
do while(n.lt.Nparticles)
    call random_number(s)
    s(1)=xmin+s(1)*(xmax-xmin)
    s(2)=xmin+s(2)*(ymax-ymin)
    if(V_x(s)<E_cut)then
        n=n+1
        x(:,n)=s(:)
    endif
enddo
!==============================================================================!
!                       Write Initial Distribution to File
!==============================================================================!
open(unit=16,file='coor_ini.dat')
do i=1,Nparticles
    write(16,*) x(:,i)
enddo
close(16)
!==============================================================================!
!                           Compute U[x_ij] 
!compute pairwise energies for the entire set of initial particle 
!no longer using single-particle potential
!==============================================================================!
do i=2,Nparticles
    do j=1,i-1
        U(i,j)=Pair_LJ_NRG(x(:,i),x(:,j))
        U(j,i)=U(i,j)
    enddo
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
    k=random_integer(1,Nparticles)
!==============================================================================!
!                   Generate coordinates for Random move trial
!random numbers generated [0,1], make it [-1,1] ==> s=2*s-1
!==============================================================================!
    call random_number(s) 
    x0=x(:,k)+mv_cutoff*(2*s-1)
!==============================================================================!
!                      Continue if V(trial) < Ecut 
!==============================================================================!
    if(V_x(x0).lt.E_cut) then
!==============================================================================!
!                   Compute Energy Change due to Trial Move
!first compute Delta V = Vnew[x_i] - V[x_i]
!Then compute delta U and add together for total change in energy
!==============================================================================!
        counter=counter+1
        U_move(k)=P_x(x0)
        Delta_E=0d0
        do j=1,Nparticles
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
        if(exp(Delta_E).ge.t1)then
           U(:,k)=U_move(:)
           U(k,:)=U_move(:)
           accept=accept+1
           x(:,k)=x0(:)
        else
!==============================================================================!
!                   Otherwise Reject Configuration
!==============================================================================!
        endif
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
open(unit=17,file='coor_fin.dat')
do i=1,Nparticles
    write(17,*) x(:,i)
enddo
!==============================================================================!
!                               output file                                    !
!==============================================================================!
write(*,*) 'Hello Universe!'
call cpu_time(t_f)
open(90,file='simulation.dat')
write(90,*) 'Nparticles ==> ', Nparticles
write(90,*) 'particle dimensionality ==> ', d
write(90,*) 'P(x) Normalization ==> ', integral_P
write(90,*) 'Number of Normalization Integrations==> ', Nsobol
write(90,*) 'Number of MMC Iterations ==> ', N_iter
write(90,*) 'c_LJ ==> ', c_LJ
write(90,*) 'Final Move Cutoff ==> ', mv_cutoff
write(*,*) 'xmin ==>', xmin
write(*,*) 'xmax ==>', xmax
write(*,*) 'ymin ==>', ymin
write(*,*) 'ymax ==>', ymax
write(90,*) 'Total Time ==> ', t_f-t_i
close(90)
end program main
