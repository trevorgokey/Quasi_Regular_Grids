PROGRAM QM_2d
  implicit none
  integer, parameter :: Nx=70,Ny=Nx
  integer :: i,j,i1,i2,j1,j2,Ndvr,IERR
  double precision, parameter :: xmin=-5.,xmax=20.,ymin=xmin,ymax=xmax,Ecut=15,pi=dacos(-1d0)
  double precision :: x(Nx),y(Ny),dx,dy,Potential,q(2),q1(2),Vpot,T1(Nx,Nx),T2(Ny,Ny)
  double precision, allocatable :: H(:,:),E(:), FV1(:),Fv2(:),z(:,:)
  open(1,file='Energies.dat')
  dx=(xmax-xmin)/(Nx-1)
  dy=(ymax-ymin)/(Ny-1)
  do i=1,Nx
     x(i)=xmin+(i-1)*dx
  enddo
  do i=1,Ny
     y(i)=ymin+(i-1)*dy
  enddo
  do i=1,Nx
     do j=1,Nx
        if(i.ne.j) then
           T1(i,j)=(-1)**(i-j)/(dx*(i-j))**2
        else
           T1(i,j)=pi**2/6/dy**2
        endif
     enddo
  enddo
  do i=1,Ny
     do j=1,Ny
        if(i.ne.j) then
           T2(i,j)=(-1)**(i-j)/(dy*(i-j))**2
        else
           T2(i,j)=pi**2/6/dy**2
        endif
     enddo
  enddo
  
  Ndvr=0
  do i=1,Nx        
     q(1)=x(i)
     do j=1,Ny
        q(2)=y(j)
        if(Potential(q)<Ecut) then
           Ndvr=Ndvr+1
        endif
     enddo
  enddo
  write(*,*) 'Nx=',Nx,' Ny=',Ny,' Ndvr=', Ndvr
  allocate(H(Ndvr,Ndvr),z(Ndvr,Ndvr),E(Ndvr),Fv1(Ndvr),Fv2(Ndvr))
  i=0
  H=0d0
  do i1=1,Nx        
     q(1)=x(i1)
     do i2=1,Ny
        q(2)=y(i2)
        Vpot=Potential(q)
        if(Vpot<Ecut) then
           i=i+1
           H(i,i)=H(i,i)+Vpot
           j=0
           do j1=1,Nx
              q1(1)=x(j1)
              do j2=1,Ny
                 q1(2)=y(j2)
                 if(Potential(q1)<Ecut) then
                    j=j+1
                    if(i1==j1) then   ! y-kinetic energy
                       if(i2.ne.j2) then
                          H(i,j)=H(i,j)+(-1)**(i2-j2)/(dy*(i2-j2))**2
                       else
                          H(i,j)=H(i,j)+pi**2/6/dy**2
                       endif
                    endif
                    if(i2==j2) then   ! x-kinetic energy
                       if(i1.ne.j1) then
                          H(i,j)=H(i,j)+(-1)**(i1-j1)/(dx*(i1-j1))**2
                       else
                          H(i,j)=H(i,j)+pi**2/6/dx**2
                       endif
                    endif
                 endif
              enddo
           enddo
        endif
     enddo
  enddo
  call RS(Ndvr,Ndvr,H,E,1,z,FV1,FV2,IERR)

  write(1,*) '#   Nx=',Nx,' Ny=',Ny,' Ecut=', Ecut,'  Ndvr=', Ndvr
  do i=1,122
     write(1,*) E(i)
  enddo
      
end PROGRAM QM_2d

function Potential(x)
!==============================================================================!
!       Discussion:
!Hard-coded Morse Potential Energy 
!==============================================================================!
  implicit none
  double precision, parameter :: D=12., a1=0.2041241, a2=0.18371169
  double precision x(2),Potential
  Potential=D*( (exp(-a1*x(1)) -1)**2 + (exp(-a2*x(2)) - 1)**2 )
  return
end function Potential
    
    
