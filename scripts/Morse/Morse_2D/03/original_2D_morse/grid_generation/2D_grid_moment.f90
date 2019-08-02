!=============================================================================80
!               2D Morse quasi-Lennard-Jones Grid Generation
!==============================================================================!
!       Discussion:
!Generate gridpoints (Particles) distributed via the 2D Morse Oscillator
!This grid is 'optimized' using a quasi-Lennard Jones Potential (12-6)
!Gridpoints have a quasi-regular distribution locally
!Gridpoints are constricted by an energy cutoff contour (E_cut)
!==============================================================================!
!       Modified:
!   19 May 2019
!       Author:
!   Shane Flynn 
!==============================================================================!
module qlj_mod
implicit none
!==============================================================================!
!                            Global Variables 
!==============================================================================!
!d              ==> Particle Dimensionality
!Npoints        ==> Number of grid points 
!c_LJ           ==> parameter for LJ
!E_cut          ==> Energy Cutoff Contour
!integral_P     ==> Normalization for P_x 
!==============================================================================!
integer::d,Npoints
double precision::c_LJ,E_cut,integral_P 
!==============================================================================!
contains
!==============================================================================!
function random_integer(Nmin,Nmax) 
!==============================================================================!
!       Discussion:
!Randomly generate an integer in the range 1-Nparticles for MMC Algorithm
!==============================================================================!
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
function V(x)
!==============================================================================!
!       Discussion:
!Potential Energy (Hard-Coded 2D Morse)
!==============================================================================!
!V              ==> evaluate V(x,y)
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!D_morse        ==> morse constant (Garaschuk 2001 for values) 
!omega          ==>(d) morse constants (Garaschuk 2001 for values)
!==============================================================================!
implicit none 
double precision::x(d),V
double precision,parameter::omega(2)=(/0.2041241,0.18371169/)
double precision,parameter::D_morse=12.
V=D_morse*sum((exp(-omega(:)*x(:))-1.)**2)
end function V
!==============================================================================!
subroutine Moments_Reg(Moment,N,xmin,xmax)
!==============================================================================!
!Compute low-level moments of the distribution to verify global accuracy
!       Discussion:
!Compute moments using a regular square grid (most accuract method for 3D Morse)
!Integrate over the square [a,b],[a,b] size determined by Moments_MMC subroutine
!int P(r)~Area_Square/N sum_n=1,N P(r_n)
!==============================================================================!
!r              ==>(d) coordinates 
!xmin           ==>(d) minimum of normalization box
!xmax           ==>(d) maximum of normalization box
!Moment         ==>(0:5) 5 Lowest Moments for the distribution
!==============================================================================!
integer::N,i1,i2,j
double precision::r(d),xmin(d),xmax(d),Moment(0:5),dummy
Moment=0d0
do i1=0,N
    do i2=0,N
        r(1)=xmin(1)+i1*(xmax(1)-xmin(1))/N
        r(2)=xmin(2)+i2*(xmax(2)-xmin(2))/N
        dummy=P(r)
        Moment(0)=Moment(0)+dummy
        do j=1,d
            Moment(j)=Moment(j)+dummy*r(j)
            Moment(d+j)=Moment(d+j)+dummy*r(j)**2
        enddo
        Moment(5)=Moment(5)+dummy*r(1)*r(2)
    enddo
enddo
dummy=1./N**d
do j=1,d
    dummy=dummy*(xmax(j)-xmin(j))
enddo
integral_P=dummy*Moment(0)
Moment(1:5)=Moment(1:5)/Moment(0)
end subroutine Moments_Reg
!==============================================================================!
subroutine Moments_MC(Moment,N_MC,xmin,xmax)
!==============================================================================!
!Compute low-level moments of the distribution to verify global accuracy
!       Discussion:
!Compute moments using Monte Carlo (not Metropolis, plain MC)
!==============================================================================!
!r              ==>(d) coordinates 
!xmin           ==>(d) minimum of normalization box
!xmax           ==>(d) maximum of normalization box
!Moment         ==>(0:5) 5 Lowest Moments for the distribution
!==============================================================================!
integer::N_MC,i,j
double precision::r(d),xmin(d),xmax(d),Moment(0:5),dummy
Moment=0d0
do i=1,N_MC
    call random_number(r)
    r(:)=xmin(:)+r(:)*(xmax(:)-xmin(:))
    dummy=P(r)
    Moment(0)=Moment(0)+dummy
    do j=1,d
        Moment(j)=Moment(j)+dummy*r(j)
        Moment(d+j)=Moment(d+j)+dummy*r(j)**2
    enddo
    Moment(5)=Moment(5)+dummy*r(1)*r(2)
enddo
dummy=1./N_MC
do i=1,d
    dummy=dummy*(xmax(i)-xmin(i))
enddo
integral_P=dummy*Moment(0)
Moment(1:5)=Moment(1:5)/Moment(0)
end subroutine Moments_MC
!==============================================================================!
subroutine Moments_MMC(Moment,N_MC,xmin,xmax)
!==============================================================================!
!Compute low-level moments of the distribution to verify global accuracy
!       Discussion:
!Compute moments using Metropolis Monte Carlo 
!This subroutine sets the box size for normalizing P (xmin,xmax)
!==============================================================================!
!N_MC           ==># of Monte Carlo Iterations
!mv_cutoff      ==>trial displacement move cutoff
!r              ==>(d) coordinates 
!r_trial        ==>(d) trial coordinates 
!s              ==>(d) trail random number for coordinates
!xmin           ==>(d) minimum of normalization box
!xmax           ==>(d) maximum of normalization box
!Moment         ==>(0:5) 5 Lowest Moments for the distribution
!==============================================================================!
integer::N_MC,i,j
double precision::Moment(0:5),dummy,r_trial(d),r(d),s(d),xmin(d),xmax(d)
double precision,parameter::mv_cutoff=0.1
Moment=0d0
r=0d0
xmin=r
xmax=r
do i=1,N_MC
!==============================================================================!
!                   Generate coordinates for Random move
!           random numbers generated (0,1), make it (-1,1) ==> s=2*s-1
!==============================================================================!
    call random_number(s) 
    r_trial=r+mv_cutoff*(2*s-1)
!==============================================================================!
!                             Test Acceptance
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
    Moment(5)=Moment(5)+r(1)*r(2)
enddo
Moment=Moment/N_MC
end subroutine Moments_MMC
!==============================================================================!
subroutine Moments(Moment,x)
!==============================================================================!
!Compute low-level moments of the distribution to verify global accuracy
!       Discussion:
!Compute moments given a set of points
!==============================================================================!
!x              ==>(d,Npoints) All coordinates 
!Moment         ==>(0:5) 5 Lowest Moments for the distribution
!==============================================================================!
integer::i,j
double precision::Moment(0:5),x(d,Npoints)
Moment=0d0
do i=1,Npoints
    do j=1,d
        Moment(j)=Moment(j)+x(j,i)
        Moment(d+j)=Moment(d+j)+x(j,i)**2
    enddo
    Moment(5)=Moment(5)+x(1,i)*x(2,i)
enddo
Moment(1:5)=Moment(1:5)/Npoints
end subroutine Moments
!==============================================================================!
function P(x)
!==============================================================================!
!       Discussion:
!Target Distribution Function (Hard-Coded)
!==============================================================================!
!P              ==>evaluate P(x)
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!Delta          ==>Delta parameter for P(x) distribution, ~10% of E_cut
!integral_P     ==>Normalization factor for P
!==============================================================================!
implicit none 
double precision::x(d),P
!double precision,parameter::Delta=0.1
if(V(x)<E_cut) then
P=(E_cut-V(x))/integral_P          
!P=(E_cut*(1+Delta)-V(x))/integral_P          
else                                              !set equal to 0 if beyond Ecut
   P=1d-20
endif
end function P
!==============================================================================!
function Pair_LJ_NRG(x1,x2)
!==============================================================================!
!       Discussion:
!quasi-LJ energy between 2 grid points
!==============================================================================!
!x1             ==>(d) ith atoms coordinates
!x2             ==>(d) jth atoms coordinates
!a/b            ==>evaluate LJ
!sigma          ==>c*sigma(P)
!Pair_LJ_NRG    ==>Energy of the i-j qLJ potential
!==============================================================================!
implicit none 
double precision::x1(d),x2(d),a,b,Pair_LJ_NRG,sigma1,sigma2
a=sum((x1(:)-x2(:))**2)
sigma1=c_LJ*(P(x1)*Npoints)**(-1./d)    
sigma2=c_LJ*(P(x2)*Npoints)**(-1./d)    
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
!==============================================================================!
!N_MC           ==># of MMC steps to distribute points
!N_MC_moments   ==># of MC moves to compute moments 
!N_reg          ==># grid points in each dimension for moments with regular grid 
!MMC_freq       ==>Interval to update MMC acceptance ratio
!Delta_E        ==>change in energy due to trial move
!mv_cutoff      ==> maximum displacement factor for trial move 
!Moment         ==>(0:5) 5 Lowest Moments for the distribution
!x              ==>(d) coordinates for ith gridpoint
!x0             ==>(d) store coordinate before trial move
!xmin           ==>(d) minimum for P(x) normalization box
!xmax           ==>(d) maximum of P(x) normalization box
!U              ==>(Npoints,Npoints) All i,j pair-wise energies 
!U_move         ==>(Npoints)  qLJ-Energy for trial movement
!==============================================================================!
implicit none
integer::N_MC,N_MC_moments,N_reg,MMC_freq,accept,counter,i,j,k,n
double precision::Delta_E,mv_cutoff,t1,time1,time2,Moment(0:5)
double precision,allocatable,dimension(:)::x0,s,U_move,xmin(:),xmax(:)
double precision,allocatable,dimension(:,:)::x,U
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
call cpu_time(time1)
read(*,*) d
read(*,*) Npoints
read(*,*) N_MC_moments  
read(*,*) N_reg     
read(*,*) N_MC     
read(*,*) MMC_freq  
read(*,*) E_cut
read(*,*) c_LJ
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(x(d,Npoints),x0(d),s(d),U(Npoints,Npoints),U_move(Npoints),xmin(d))
allocate(xmax(d))
write(*,*) 'Test 0; Successfully Read Input File'
!==============================================================================!
!                   Run Monte Carlo to Normalize P/Get Moments
!==============================================================================!
integral_P=1d0                       !set equal to 1 so you can initially call P
call Moments_MMC(Moment,N_MC_moments,xmin,xmax)
write(*,*) 'Moments from Metropolis MC usig N_MC steps=',N_MC_moments
write(*,*) Moment(1:5)
do i=1,d
    write(*,*) i,' xmin xmax ==>', xmin(i),xmax(i)
enddo
xmin=xmin-(xmax-xmin)*0.01
xmax=xmax+(xmax-xmin)*0.01
call Moments_Reg(Moment,N_reg,xmin,xmax)
write(*,*) 'Moments using regular integration:'
write(*,*) Moment(1:5)
write(*,*) 'Normalization of P =', Integral_P
write(*,*) 'Test 1; Successfully normalized P(x)'
!==============================================================================!
!                       Generate Initial Distribution
!               Initally accept any point where Potential<Ecut
!==============================================================================!
i=1
do while(i.le.Npoints)   
    call random_number(s)
    s(:)=xmin(:)+s(1)*(xmax(:)-xmin(:))
    if(V(s)<E_cut)then
        x(:,i)=s(:)
        i=i+1            
    endif
enddo
!==============================================================================!
!                          Write Initial Coordinates
!==============================================================================!
open(unit=17,file='coor_ini.dat')
do i=1,Npoints
    write(17,*) x(:,i)
enddo
close(17)
call Moments(Moment,x)
write(*,*) 'Moments from the initial distribution, # Gridpoints ==>',Npoints
write(*,*) Moment(1:5)
write(*,*) 'Test 2; Successfully generated initial grid' 
!==============================================================================!
!                           Compute U[x_ij] 
!             Pairwise energies for all the initial Grid Points
!==============================================================================!
do i=2,Npoints
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
do n=1,N_MC
!==============================================================================!
!                           Select Atom to Move
!==============================================================================!
    k=random_integer(1,Npoints)
!==============================================================================!
!                           Generate trial move 
!           random numbers generated (0,1), make it (-1,1) ==> s=2*s-1
!==============================================================================!
    call random_number(s) 
    x0=x(:,k)+mv_cutoff*(2*s-1)
!==============================================================================!
!                   Only consider point if V(trial)<Ecut 
!                 Compute Energy Change due to Trial Move
!==============================================================================!
    if(V(x0).lt.E_cut) then
        counter=counter+1
        U_move(k)=P(x0)
        Delta_E=0d0
        do j=1,Npoints
            if(j.ne.k) then
                U_move(j)=Pair_LJ_NRG(x(:,j),x0)
                Delta_E=Delta_E+U(j,k)-U_move(j)
            endif
        enddo
!==============================================================================!
!               Test to see if you should accept higher energy move
!           if delta e > 0 than always larger than 1, always accept
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
!       acceptance ~50%, adjust trial movement displacement accordingly
!==============================================================================!
        if(mod(n,MMC_freq)==0)then
            write(*,*) 'MMC Iteration', n
            if(dble(accept)/counter<0.5)then 
                mv_cutoff=mv_cutoff*0.9
            else
                mv_cutoff=mv_cutoff*1.1
            endif
        accept=0
        counter=0
        call Moments(Moment,x)
        write(*,*) 'MMC Moments:'
        write(*,*) Moment(1:5)
        endif
enddo
!==============================================================================!
!                           Write Optimized Grid 
!==============================================================================!
open(unit=18,file='grid.dat')
do i=1,Npoints
    write(18,*) x(:,i)
enddo
close(18)
write(*,*) 'Test 3; Successfully Generated Quasi-Regular Grid'
!==============================================================================!
!                               Output File                                    !
!==============================================================================!
call cpu_time(time2)
open(90,file='simulation.dat')
write(90,*) 'particle dimensionality ==> ', d
write(90,*) 'Npoints ==> ', Npoints
write(90,*) 'P Normalization Iterations (N_MC) ==> ', N_MC_moments
write(90,*) 'P Normalization factor ==> ', integral_P
write(90,*) 'Number of MMC Iterations for Grid ==> ', N_MC
write(90,*) 'c_LJ ==> ', c_LJ
write(90,*) 'E_cut ==> ', E_cut
do i=1,d
    write(90,*) i,' xmin xmax ==>', xmin(i),xmax(i)
enddo
write(90,*) 'Total Time ==> ', time2-time1
close(90)
write(*,*) 'Test 4; Hello Universe!'
end program main
