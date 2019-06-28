function integral_P(NMC)
  integer :: NMC,n
  double precision :: integral_P,xmin,xmax,ymin,ymax,x(2)
  xmin=-4
  xmax=20
  ymin=-4
  ymax=20
  integral_P=0d0
  do n=1,NMC
     call random_number(x)
     x(1)=xmin+x(1)*(xmax-xmin)
     x(2)=xmin+x(2)*(ymax-ymin)
     integral_P=integral_P+P(x)
  enddo
  integral_P=integral_P*(xmax-xmin)*(ymax-ymin)/NMC
end function integral_P
     
     
