!=============================================================================80
!                       2D Henon-Heiles QRG Generation
!==============================================================================!
!       Discussion:
!Generate gridpoints distributed via the 2D Henon-Helies Potential
!This QRG is using a quasi-Lennard Jones Potential
!An energy cutoff contour is defined at 13.3333 (bound states of the potential)
!Force minimization of the gridpoints, accepting any trial moves if it decreases
!the energy.
!Compute the lower moments of the distribution to test accuracy
!Moments are computed in a few ways for relative comparisons.
!==============================================================================!
!       Modified:
!   24 June 2019
!       Author:
!   Shane Flynn
!==============================================================================!
module qlj_mod
implicit none
!==============================================================================!
!                            Global Variables
!==============================================================================!
!d              ==>Particle Dimensionality (hard-coded Potential 2D only)
!Npoints        ==>Number of grid points
!c_LJ           ==>parameter for qLJ (should be on the order of 1)
!E_cut          ==>Energy Cutoff Contour
!integral_P     ==>Normalization for P_x (computed with the moments)
!==============================================================================!
integer,parameter::d=2
integer::Npoints
double precision::c_LJ,E_cut,integral_P
!==============================================================================!
contains
!==============================================================================!
function random_integer(Nmin,Nmax)
!==============================================================================!
!       Discussion:
!Pseudo-Random generation of an integer in the range Nmin:Nmax
!==============================================================================!
!Nmin/max       ==>Range to generate integers from
!random_integer ==>Integer returned
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
!Potential Energy (Hard-Coded 2-Dimensional Henon-Heiles)
!==============================================================================!
!V              ==>evaluate V(x,y)
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!lambda         ==>HH Potential Parameter: Hamilton,Light 1986 paper for values
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
!Target Distribution Function. The form comes from a simple semi-classical
!argument proposed by Bill Poirier
!Found. Phys. 30, 1191 (2000) [specifically Eq.(12)]
!==============================================================================!
!P              ==>evaluate P(x), the target distribution function
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!integral_P     ==>Normalization constant for P(x),(see moments subroutines)
!==============================================================================!
implicit none
double precision::x(d),P
if(V(x)<E_cut) then
  P=(E_cut-V(x))/integral_P
else                                              !set equal to 0 if beyond Ecut
  P=1d-20
endif
end function P
!==============================================================================!
function Pair_qLJ_NRG(x1,x2)
!==============================================================================!
!       Discussion:
!Pseudo-Pairwise Potential for minimizing distances between grid (quasi-LJ)
!==============================================================================!
!x1             ==>(d) ith atoms coordinates
!x2             ==>(d) jth atoms coordinates
!a/b            ==>evaluate qLJ
!sigma          ==>inverse gaussian widths
!Pair_LJ_NRG    ==>Energy of the i-j q-LJ potential
!==============================================================================!
implicit none
double precision::x1(d),x2(d),a,b,Pair_qLJ_NRG,sigma1,sigma2
a=sum((x1(:)-x2(:))**2)
sigma1=c_LJ*(P(x1)*Npoints)**(-1./d)
sigma2=c_LJ*(P(x2)*Npoints)**(-1./d)
b=(sigma2**2/a)**3
a=(sigma1**2/a)**3
Pair_qLJ_NRG=a**2-a+b**2-b
end function Pair_qLJ_NRG
!==============================================================================!
subroutine Moments_Reg(Moment,i,xmin,xmax)
!==============================================================================!
!       Discussion:
!Compute moments of the distribution (to verify global accuracy)
!Compute moments using a regular square grid (most accurate method for 3D Morse)
!Integrate over the square [a,b],[a,b] size determined by Moments_MMC subroutine
!int P(r)~Area_Square/N sum_n=1,N P(r_n)
!==============================================================================!
!r              ==>(d) coordinates
!xmin           ==>(d) minimum of normalization square grid
!xmax           ==>(d) maximum of normalization square grid
!Moment         ==>(0:5) 5 Lowest Moments for the distribution
!==============================================================================!
integer::i,i1,i2,j
double precision::r(d),xmin(d),xmax(d),Moment(0:5),dummy
Moment=0d0
do i1=0,i
  do i2=0,i
    r(1)=xmin(1)+i1*(xmax(1)-xmin(1))/i
    r(2)=xmin(2)+i2*(xmax(2)-xmin(2))/i
    dummy=P(r)
    Moment(0)=Moment(0)+dummy
    do j=1,d
      Moment(j)=Moment(j)+dummy*r(j)
      Moment(d+j)=Moment(d+j)+dummy*r(j)**2
    enddo
    Moment(5)=Moment(5)+dummy*r(1)*r(2)
  enddo
enddo
dummy=1./i**d
do j=1,d
  dummy=dummy*(xmax(j)-xmin(j))
enddo
integral_P=dummy*Moment(0)
Moment(1:5)=Moment(1:5)/Moment(0)
end subroutine Moments_Reg
!==============================================================================!
subroutine Moments_MC(Moment,N_MC,xmin,xmax)
!==============================================================================!
!       Discussion:
!Compute moments of the distribution (to verify global accuracy)
!Compute moments using Monte Carlo (Monte Carlo, not Metropolis)
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
!       Discussion:
!Compute moments of the distribution (to verify global accuracy)
!Compute moments using Metropolis Monte Carlo
!This subroutine determines the box size for normalizing P; (xmin,xmax)
!==============================================================================!
!N_MC           ==>Number of Monte Carlo Iterations
!mv_cutoff      ==>trial move-displacement cutoff
!r              ==>(d) coordinates
!r_trial        ==>(d) trial coordinates
!s              ==>(d) Pseudo-randomly move coordinates for trial move
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
      if(xmin(j)>r(j)) then
        xmin(j)=r(j)
      elseif (xmax(j)<r(j)) then
        xmax(j)=r(j)
      endif
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
!       Discussion:
!Compute moments of the distribution (to verify global accuracy)
!Compute moments directly from a set of points
!==============================================================================!
!x              ==>(d,Npoints) All Coordinates
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
end module qlj_mod
!==============================================================================!
!==============================================================================!
!==============================================================================!
program main
use qlj_mod
!==============================================================================!
!==============================================================================!
!==============================================================================!
!N_MC           ==># of MC steps to distribute points
!N_MC_moments   ==># of MC moves to compute moments
!N_reg          ==># grid points (per dimension) for moments with regular grid
!MMC_freq       ==>Interval to update MC acceptance ratio
!accept         ==># of acceptances for adjusting trial move parameter
!Delta_E        ==>change in energy due to trial move
!mv_cutoff      ==>maximum displacement for trial move
!Moment         ==>(0:5) 5 Lowest Moments of the distribution
!x0             ==>(d) store coordinate before trial move
!U_move         ==>(Npoints) qLJ-Energy for trial movement
!xmin           ==>(d) minimum of P(x) normalization box
!xmax           ==>(d) maximum of P(x) normalization box
!x              ==>(d,Npoints) All coordinates
!U              ==>(Npoints,Npoints) All i,j pair-wise energies
!==============================================================================!
implicit none
integer::N_MC,N_MC_moments,N_reg,MMC_freq,accept,counter,i,j,k,plt_count
double precision::Delta_E,deltae1,mv_cutoff,time1,time2,Moment(0:5)
double precision,allocatable,dimension(:)::x0,s,U_move,xmin(:),xmax(:)
double precision,allocatable,dimension(:,:)::x,U
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
call cpu_time(time1)
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
integral_P=1d0                      !set equal to 1 so you can initially call P
call Moments_MMC(Moment,N_MC_moments,xmin,xmax)
write(*,*) 'Moments from Metropolis MC usig N_MC steps=',N_MC_moments
write(*,*) Moment(1:5)
do i=1,d
  write(*,*) i,' xmin, xmax ==>', xmin(i),xmax(i)
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
  s(:)=xmin(:)+s(:)*(xmax(:)-xmin(:))
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
write(*,*) 'Moments from the initial distribution, # Gridpoints ==>', Npoints
write(*,*) Moment(1:5)
write(*,*) 'Test 2; Successfully generated initial grid'
!==============================================================================!
!                           Compute U[x_ij]
!             Pairwise energies for all the initial Grid Points
!==============================================================================!
do i=2,Npoints
  do j=1,i-1
    U(i,j)=Pair_qLJ_NRG(x(:,i),x(:,j))
    U(j,i)=U(i,j)
  enddo
enddo
!==============================================================================!
!                       Begin MMC to Optimize GridPoints
!==============================================================================!
accept=0
counter=0
mv_cutoff=0.01
deltae1=0
plt_count=0
open(unit=18,file='mv_cut.dat')
open(unit=19,file='delE.dat')
do i=1,N_MC
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
        U_move(j)=Pair_qLJ_NRG(x(:,j),x0)
        Delta_E=Delta_E+U(j,k)-U_move(j)
      endif
    enddo
!==============================================================================!
!          Accept any Trial that decreases the energy (force minimizaton)
!==============================================================================!
    if(Delta_E.ge.0d0)then
      U(:,k)=U_move(:)
      U(k,:)=U_move(:)
      accept=accept+1
      x(:,k)=x0(:)
      deltae1=deltae1+Delta_E
    endif
  endif
!==============================================================================!
!                           Update Cutoff Paramater
!       acceptance ~50%, adjust trial movement displacement accordingly
!==============================================================================!
  if(mod(i,MMC_freq)==0)then
    write(*,*) 'MMC Iteration', i
    if(dble(accept)/counter<0.3)then
      mv_cutoff=mv_cutoff*0.9
    else
      mv_cutoff=mv_cutoff*1.1
    endif
    accept=0
    counter=0
    call Moments(Moment,x)
    write(*,*) 'MMC Moments:'
    write(*,*) Moment(1:5)
    write(*,*) 'mv cutoff', mv_cutoff
    write(*,*) 'Deltae1==>', Deltae1
    plt_count=plt_count+1
    write(18,*) plt_count, mv_cutoff
    write(19,*) plt_count, Deltae1
    deltae1=0
  endif
enddo
close(18)
close(19)
!==============================================================================!
!                           Write Optimized Grid
!==============================================================================!
open(unit=20,file='QRG.dat')
do i=1,Npoints
  write(20,*) x(:,i)
enddo
close(20)
write(*,*) 'Test 3; Successfully Generated Quasi-Regular Grid'
!==============================================================================!
!                               Output File                                    !
!==============================================================================!
call cpu_time(time2)
open(99,file='simulation.dat')
write(99,*) 'Grid dimensionality ==> ', d
write(99,*) 'Npoints ==> ', Npoints
write(99,*) 'P Normalization Iterations (N_MC) ==> ', N_MC_moments
write(99,*) 'P Normalization factor ==> ', integral_P
write(99,*) 'Number of MMC Iterations for Grid ==> ', N_MC
write(99,*) '# gridpoints regular moments', N_reg
write(99,*) '# MC moves to distribute points', N_MC
write(99,*) 'c_LJ ==> ', c_LJ
write(99,*) 'E_cut ==> ', E_cut
write(99,*) 'final mv_cutoff==> ', mv_cutoff
write(99,*) 'MMC frequency', MMC_freq
do i=1,d
  write(99,*) i,' xmin xmax ==>', xmin(i),xmax(i)
enddo
write(99,*) 'Total Time ==> ', time2-time1
close(99)
write(*,*) 'Hello Universe!'
end program main
