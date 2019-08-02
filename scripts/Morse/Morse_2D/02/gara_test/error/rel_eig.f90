program main
implicit none
!read in relative errors and plot an average over every 5 eigenvalues
integer::i,neig
double precision,allocatable::x(:)

neig=120
allocate(x(neig+1))
open(15,File="relative.dat")
read(15,*)x
close(15)

!write(*,*) 'x'
!write(*,*) x

do i=1,neig
write(*,*) i, x(i+1)
enddo

end program main
