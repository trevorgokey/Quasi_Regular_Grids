
module moments_mod
implicit none
use util
contains
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
end module moments_mod
