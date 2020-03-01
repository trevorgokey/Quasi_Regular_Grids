
!==============================================================================!
!==============================================================================!
program main
use qlj_mod
use omp_lib
use gridgen
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

call gridgen_run()
write(*,*) 'Hello Universe!'
end program main
