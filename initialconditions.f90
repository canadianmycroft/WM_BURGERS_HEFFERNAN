program initialconditions
implicit none
!this program provides the initial conditions for various runs of the burgers equation program
real :: f_x,l,g_x
integer :: x,n
open(1,file='initialconditions')

print *,'Enter an integer step size'
read *, n

do x=1,n
	f_x=sin((2*3.14159265*x)/n)
	write(1,*) f_x
end do
end program initialconditions