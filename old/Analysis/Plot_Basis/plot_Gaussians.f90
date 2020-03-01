program plot_Gaussians
  implicit none
  integer, parameter :: N=200
  integer i,j
  double precision, parameter :: pi=dacos(-1d0)
  double precision :: sigma, x, y
  open(1,file='200.dat')
  open(2,file='200gaussians.dat')
  do i=1,N
     read(1,*) x,y,sigma
     sigma=1/sqrt(sigma)
     do j=1,101
        write(2,*) x+sigma*cos(2.*pi*j/100.), y+sigma*sin(2.*pi*j/100.)
     enddo
     write(2,*)
  enddo
end program plot_Gaussians
