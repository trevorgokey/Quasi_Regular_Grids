!=============================================================================80
!                    3D Morse Quasi-Reglar Grid Generation
!==============================================================================!
!Generate gridpoints distributed via the 3D Morse Oscillator
!This Quasi-Regular gird is optimized using a quasi-Lennard Jones Potential
!Force minimization, accept any trial moves that reduces the system's energy
!Compute low-level moments to check distribiution accuracy
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
!d              ==>Particle Dimensionality
!Npoints        ==>Number of grid points
!c_LJ           ==>parameter for LJ
!E_cut          ==>Energy Cutoff Contour
!integral_P     ==>Normalization for P_x
!==============================================================================!
integer::d,Npoints,Natoms
double precision::c_LJ,E_cut,integral_P
double precision,allocatable,dimension(:,:)::R
!==============================================================================!
contains
!==============================================================================!

function load_density( fnm)
    implicit none
    integer::Natoms,i,load_density
    character(len = *) :: fnm
    open( 15, file=fnm)
    read( 15, *) Natoms
    allocate (R(d, Natoms))
    read( 15, *)
    do i = 1, Natoms
        read(15, *) R(:,i)
    end do
    close( 15)
    load_density = Natoms
end function load_density

function random_integer(Nmin,Nmax, a) 
!==============================================================================!
!Randomly generate an integer in the range 1-Nparticles
!==============================================================================!
!Nmin           ==>minimum index value (1)
!Nmax           ==>maximum index value (Nparticles)
!random_integer ==>integer returned
!a              ==>random number (0,1)
!==============================================================================!
implicit none
integer::Nmin,Nmax,random_integer
double precision::a
!call random_number(a)
random_integer=floor(a*(Nmax-Nmin+1))+Nmin
end function random_integer
!==============================================================================!
function V(x)
!==============================================================================!
!Potential Energy (Hard-Coded 3D Morse)
!==============================================================================!
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!V              ==>evaluate V(x)
!D_morse        ==>Parameter for Morse Potential
!omega          ==>(d) Parameter for Morse Potential
!==============================================================================!
implicit none 
integer::i
double precision::x(d),V
double precision,parameter::omega(3)=(/0.2041241,0.18371169,0.16329928/)
double precision,parameter::D_morse=12.
double precision ::  rho=0., sigma=30., dist

!open(unit=40,file='grid.xyz', status="REPLACE")
!close(40)

rho = 0.0
do i = 1, Natoms
dist = sum((x(:) - R(:,i))**2)
!if ( dist < sigma*3 ) then
    rho = rho - exp(-dist/sigma) 
!end if
end do
V = rho
end function V
!==============================================================================!
function P(x)
!==============================================================================!
!Target Distribution Function
!==============================================================================!
!P              ==>evaluate P(x)
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!==============================================================================!
implicit none 
double precision::x(d),P
if(V(x)<E_cut) then
   P=(E_cut-V(x))/integral_P          
else                                             !set equal to 0 if beyond Ecut
   P=1d-20
end if
end function P
!==============================================================================!
function Pair_LJ_NRG(x1,x2)
!==============================================================================!
!quasi-Lennard Jones pairwise energy between grid points
!==============================================================================!
!x1             ==>(d) ith atoms coordinates
!x2             ==>(d) jth atoms coordinates
!a/b            ==>evaluate LJ
!sigma          ==>c*sigma(P)
!Pair_LJ_NRG    ==>Energy of the i-j q-LJ potential
!==============================================================================!
implicit none 
double precision::x1(d),x2(d),a,b,Pair_LJ_NRG,sigma1,sigma2
a=sum((x1(:)-x2(:))**2)
Pair_LJ_NRG=0.0
!if ( a < 20**2 ) then
    sigma1=c_LJ*(P(x1)*Npoints)**(-1./d)    
    sigma2=c_LJ*(P(x2)*Npoints)**(-1./d)    
    b=(sigma2**2/a)**3
    a=(sigma1**2/a)**3
    Pair_LJ_NRG=a**2-a+b**2-b
!end if
Pair_LJ_NRG = Pair_LJ_NRG
end function Pair_LJ_NRG
!==============================================================================!
subroutine Moments_Reg(Moment,i,xmin,xmax)
!==============================================================================!
!Compute lower moments of the distribution to verify global accuracy
!Compute moments using a regular square grid (most accurate method for 3D Morse)
!Integrate over the square [a,b],[a,b] size determined by Moments_MMC subroutine
!int P(r)~Area_Square/N sum_n=1,N P(r_n)
!==============================================================================!
!r              ==>(d) coordinates
!xmin           ==>(d) minimum of normalization box
!xmax           ==>(d) maximum of normalization box
!Moment         ==>(0:5) 5 Lowest Moments for the distribution
!==============================================================================!
integer::i,i1,i2,i3,j
double precision::r(d),xmin(d),xmax(d),Moment(0:9),dummy
Moment=0d0
do i1=0,i
    do i2=0,i
        do i3=0,i
            r(1)=xmin(1)+i1*(xmax(1)-xmin(1))/i
            r(2)=xmin(2)+i2*(xmax(2)-xmin(2))/i
            r(3)=xmin(3)+i3*(xmax(3)-xmin(3))/i
            !dummy=P(r)
            !Moment(0)=Moment(0)+dummy
            !do j=1,d
            !    Moment(j)=Moment(j)+dummy*r(j)
            !    Moment(d+j)=Moment(d+j)+dummy*r(j)**2
            !enddo
            !Moment(7)=Moment(7)+dummy*r(1)*r(2)
            !Moment(8)=Moment(8)+dummy*r(2)*r(3)
            !Moment(9)=Moment(9)+dummy*r(3)*r(1)
        enddo
    enddo
enddo
!dummy=1./i**3
!do j=1,d
!   dummy=dummy*(xmax(j)-xmin(j))
!enddo
!integral_P=dummy*Moment(0)
!Moment(1:9)=Moment(1:9)/Moment(0)
end subroutine Moments_Reg
!==============================================================================!
subroutine Moments_MC(Moment,N_MC,xmin,xmax)
!==============================================================================!
!==============================================================================!
!Compute low-level moments of the distribution to verify global accuracy
!Compute moments using Monte Carlo (not Metropolis, just plain MC)
!==============================================================================!
!r              ==>(d) coordinates
!xmin           ==>(d) minimum of normalization box
!xmax           ==>(d) maximum of normalization box
!Moment         ==>(0:9) 5 Lowest Moments for the distribution
!==============================================================================!
integer::N_MC,i,j
double precision::r(d),xmin(d),xmax(d),Moment(0:9),dummy
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
   Moment(7)=Moment(7)+dummy*r(1)*r(2)
   Moment(8)=Moment(8)+dummy*r(2)*r(3)
   Moment(9)=Moment(9)+dummy*r(3)*r(1)
enddo
dummy=1./N_MC
do i=1,d
   dummy=dummy*(xmax(i)-xmin(i))
enddo
integral_P=dummy*Moment(0)
Moment(1:9)=Moment(1:9)/Moment(0)
end subroutine Moments_MC
!==============================================================================!
subroutine Moments_MMC(Moment,N_MC,xmin,xmax)
!==============================================================================!
!Compute low-level moments of the distribution to verify global accuracy
!Compute moments using Metropolis Monte Carlo 
!This subroutine determins the box size for normalizing P; (xmin,xmax)
!==============================================================================!
!N_MC           ==>Number of Monte Carlo Iterations
!mv_cutoff      ==>trial displacement move cutoff
!r              ==>(d) coordinates 
!r_trial        ==>(d) trial coordinates 
!s              ==>(d) trail displacement; random number for coordinates
!xmin           ==>(d) minimum of normalization box
!xmax           ==>(d) maximum of normalization box
!Moment         ==>(0:9) 9 Lowest Moments for the distribution
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
    !do j=1,d
    !    Moment(j)=Moment(j)+r(j)
    !    Moment(d+j)=Moment(d+j)+r(j)**2
    !enddo
    !Moment(7)=Moment(7)+r(1)*r(2)
    !Moment(8)=Moment(8)+r(2)*r(3)
    !Moment(9)=Moment(9)+r(3)*r(1)
enddo
!Moment=Moment/N_MC
end subroutine Moments_MMC
!==============================================================================!
subroutine Moments(Moment,x)
!==============================================================================!
!Compute low-level moments of the distribution to verify global accuracy
!Compute moments given a set of points.
!==============================================================================!
!x              ==>(d,Npoints) All Coordinates
!Moment         ==>(0:9) 9 Lowest Moments for the distribution
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
!==============================================================================!
end module qlj_mod
!==============================================================================!
!==============================================================================!
program main
use qlj_mod
use omp_lib
!==============================================================================!
!==============================================================================!
!N_MC           ==># of MMC steps to distribute points
!N_MC_moments   ==># of MC moves to compute moments
!N_reg          ==># grid points in each dimension for moments with regular grid
!MMC_freq       ==>Interval to update MMC acceptance ratio
!accept         ==># of acceptances for adjusting trial move acceptance
!Delta_E        ==>change in energy due to trial move
!mv_cutoff      ==>maximum displacement factor for trial move
!Moment         ==>(0:5) 5 Lowest Moments for the distribution
!x0             ==>(d) store coordinate before
!U_move         ==>(Npoints)  qLJ-Energy for trial movement
!xmin           ==>(d) minimum for P(x) normalization box
!xmax           ==>(d) maximum of P(x) normalization box
!x              ==>(d,Npoints) All coordinates
!U              ==>(Npoints,Npoints) All i,j pair-wise energies
!==============================================================================!
implicit none
integer::N_MC,N_MC_moments,N_reg,MMC_freq,accept,counter,i,j,k,plt_count,readin
integer::limit
double precision::Delta_E,deltae1,mv_cutoff,time1,time2,Moment(0:9),diff
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
Natoms = load_density('input.xyz')
write(*,*) "Natoms is", Natoms
write(*,*) 'Test 0; Successfully Read Input File and Density'
!==============================================================================!
!                 Run Monte Carlo to Normalize P/Get Moments!
!==============================================================================!
integral_P=1d0                       !set equal to 1 so you can initially call P
!call Moments_MMC(Moment,N_MC_moments,xmin,xmax)
!write(*,*) 'Moments from Metropolis MC using N_MC=',N_MC_moments
!write(*,*) Moment(1:9)
xmin=xmin-(xmax-xmin)*0.01
xmax=xmax+(xmax-xmin)*0.01
!==============================================================================!
!                       Generate Initial Distribution
!                Initally accept any point where Potential<Ecut
!==============================================================================!
!i=1
!do while(i.le.Npoints)   
!    call random_number(s)
!    write(*,*) "dist step", i, Npoints, s
!    s(:)=xmin(:)+s(:)*(xmax(:)-xmin(:))
!    write(*,*) "dist step", i, Npoints, s
!    if(V(s)<E_cut)then
!        x(:,i)=s(:)
!        i=i+1            
!    endif
!enddo
open( 15, file='input.xyz')
read( 15, *) Natoms
read( 15, *)
write(*,*) "Natoms is after read2", Natoms
k = 0
do i = 1, Npoints
k = k + 1
if (k <= Natoms) then
    read(15, *) x(:,i)
    readin = readin + 1
    write(*,*) "Read i=", i, "k=",k,Natoms, Npoints,"ik=",(i-1)*Natoms/Npoints, x(:,i)
    do j=1 , Natoms/Npoints - 1
        k = k+1
        if (k < Natoms) then
            write(*,*) "reading", i, j, k, "steps", nint(dble(Natoms)/dble(Npoints)) - 1
            read(15, *)
        else if ( k == Natoms) then
            write(*,*) "keeping last", i, j, k, "steps", nint(dble(Natoms)/dble(Npoints)) - 1
            read(15, *) x(:,i+1)
            readin = readin + 1
        else
            write(*,*) "skipping", i, j, k
        endif 
    enddo
endif
end do
close( 15)

xmin(1) = minval(x(:,1))*1.2
xmin(2) = minval(x(:,2))*1.2
xmin(3) = minval(x(:,3))*1.2
xmax(1) = maxval(x(:,1))*1.2
xmax(2) = maxval(x(:,2))*1.2
xmax(3) = maxval(x(:,3))*1.2
do i=1,d
    write(*,*) i,' xmin xmax ==>', xmin(i),xmax(i)
enddo
Integral_P = 1.0
do i=1,d
    Integral_P = Integral_P * (abs(xmax(i) - xmin(i)))
enddo
!call Moments_Reg(Moment,N_reg,xmin,xmax)
!write(*,*) 'Moments using regular integration:'
!write(*,*) Moment(1:9)
write(*,*) 'Integral of P =', Integral_P
write(*,*) 'Test 1; Successfully normalized P(x)'
open(unit=18,file='coor_ini.xyz')
write(18,*) Npoints
write(18,*) 
do i=1,Npoints
    write(18,*) "C    ", x(:,i)
enddo
close(18)

open(unit=17,file='coor_ini.dat')
do i=1,Npoints
    write(17,*) x(:,i)
enddo
close(17)
write(*,*) 'Test 2; Successfully Generated Initial GridPoints'
!call Moments(Moment,x)
!write(*,*) 'Moments from the initial distribution, Npoints= ',Npoints
!write(*,*) Moment(1:9)
!==============================================================================!
!                           Compute U[x_ij] 
!                   Pairwise energies for initial Grid
!==============================================================================!
!$OMP parallel default(none) shared(Npoints,U,x,R,Natoms) private(i,j)
!$OMP do 
do i=2,Npoints
    do j=1,i-1
        write(*,*) i,j,Npoints,i-1
        U(i,j)=Pair_LJ_NRG(x(:,i),x(:,j))
        U(j,i)=U(i,j)
    enddo
enddo
!$OMP end do
!$OMP end parallel
!==============================================================================!
!                       Begin MMC to Optimize GridPoints
!==============================================================================!
accept=0
counter=0
mv_cutoff=5
deltae1=0
plt_count=0
open(unit=18,file='mv_cut.dat')
open(unit=19,file='delE.dat')
do i=1,N_MC
write(*,*) i,N_MC
do k=1, Npoints
!==============================================================================!
!                       Randomly Select Atom to Move
!==============================================================================!
    call random_number(s) 
    !k=random_integer(1,Npoints, s(1))
!==============================================================================!
!                           Generate trial move
!        random numbers generated (0,1), make it (-1,1) ==> s=2*s-1
!==============================================================================!
    x0=x(:,k)+mv_cutoff*(2*s-1)
!==============================================================================!
!                   Only consider point if V(trial) < Ecut 
!                   Compute Energy Change due to Trial Move
!==============================================================================!
    if(V(x0).lt.E_cut) then
        counter=counter+1
        Delta_E=0d0
        U_move(k)=P(x0)
!$OMP parallel default(none) shared(Npoints,U,x,k,U_move,R,Natoms,x0) private(j,diff) reduction(+:Delta_E)
!$OMP do 
        do j=1,Npoints
            !write(*,*) "k in parallel is", k, "j is", j
            if(j.ne.k) then
                U_move(j)=Pair_LJ_NRG(x(:,j),x0)
                Delta_E = Delta_E+ ( U(j,k) - U_move(j))
                !write(*,*) j,k,"Ujk=", U(j,k), "trial U_move was", U_move(j), "leading to ", U(j,k)-U_move(j)
            endif
        enddo
!$OMP end do
!$OMP end parallel
        !write(*,*) "Delta_E from reduce is", Delta_E
!==============================================================================!
!                  Accept any trial that decreases the energy 
!==============================================================================!
        !write(*,*) "new is", U_move(k), "old was", U(k,k), "P(x) is", P(x0)
        if(Delta_E>0d0)then
            U(:,k)=U_move(:)
            U(k,:)=U_move(:)
            accept=accept+1
            x(:,k)=x0(:)
            deltae1=deltae1+Delta_E
            write(*,*) k, Delta_E, "Accepted! Displacement was", sqrt(sum((mv_cutoff*(2*s-1))**2)), "V(x0)", V(x0) 
        endif
    endif
!==============================================================================!
!                           Update Cutoff Paramater
!           acceptance rate ~50%, adjust trial displacement accordingly
!==============================================================================!
    enddo
    if(mod(i,MMC_freq)==0)then
        write(*,*) 'MMC Iteration', i
        if(dble(accept)/counter<0.5)then 
            mv_cutoff=mv_cutoff*0.9
        else
            mv_cutoff=mv_cutoff*1.1
        endif
        if (mv_cutoff < .01) then
            mv_cutoff = .01
        endif
        accept=0
        counter=0
        !call Moments(Moment,x)
        !write(*,*) 'Moments:'
        !write(*,*) Moment(1:9)
        write(*,*) 'mv cutoff', mv_cutoff
        write(*,*) 'Deltae1==>', Deltae1
        plt_count=plt_count+1
        write(18,*) plt_count, mv_cutoff
        write(19,*) plt_count, Deltae1
        deltae1=0
    endif
    open(unit=20,file='grid.xyz', access='APPEND')

        write(20,*) Npoints
        write(20,*)
        do j=1,Npoints
            write(20,*) "C     ", x(:,j)
        enddo
    close(20)
    !write(*,*) "wrote xyz"
    !open(unit=21,file='grid.dat')
    !    do j=1,Npoints
    !        write(21,*) x(:,j)
    !    enddo
    !close(21)
enddo
close(18)
close(19)
!==============================================================================!
!                         Write Quasi-Regular Grid
!==============================================================================!
open(unit=20,file='grid.dat')
do i=1,Npoints
    write(20,*) x(:,i)
enddo
close(20)
write(*,*) 'Test 3; Successfully Generated Quasi-Regular Grid'
!==============================================================================!
!                               output file                                    !
!==============================================================================!
call cpu_time(time2)
open(99,file='simulation.dat')
write(99,*) 'particle dimensionality ==> ', d
write(99,*) 'Npoints ==> ', Npoints
write(99,*) 'P Normalization Iterations (N_MC) ==> ', N_MC_moments
write(99,*) 'P(x) Normalization ==> ', integral_P
write(99,*) 'c_LJ ==> ', c_LJ
write(99,*) 'move cutoff==> ', mv_cutoff
write(99,*) 'E_cut ==> ', E_cut
write(99,*) 'N_MC_moments', N_MC_moments  
write(99,*) '# gridpoints regular moments', N_reg     
write(99,*) '# MC moves to distribute points', N_MC     
write(99,*) 'MMC frequency', MMC_freq  
do i=1,d
 write(99,*) i,' xmin xmax ==>', xmin(i),xmax(i)
enddo
write(99,*) 'Total Time ==> ', time2-time1
close(99)
write(*,*) 'Hello Universe!'
end program main
