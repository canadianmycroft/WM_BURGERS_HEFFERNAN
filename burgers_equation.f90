module globals
	implicit none
	real,allocatable :: u(:),x(:),h(:),u_old(:),numerical_flux(:)
	real :: a,b,D,u_l,u_r,u_star
	integer :: n !number of gridpoints
	character(len=1) :: what_burgers_equation,initial_data,what_continuous
	integer :: what_method
	integer :: output=81,input=80
	integer :: margin
	integer :: n_intv
	real :: n_intv_i
	integer :: i,j,k,quarter,half
	real,parameter :: pi=3.1415927
	real :: courant_number,delta_t,delta_x,time,t_max
	logical :: l_per
	logical :: first_time=.true.
	real :: f_plus,f_minus
	real :: sigma,a_roof,u1,u2,minmod
	
	contains
		function f(u)
			real :: f,u
			f=0.5*u**2
		end function
		function f_der(u)
			real :: f_der,u
			f_der=u
		end function

end module globals

program burgers_equation
	!this program computes numerical solutions to the inviscid and viscid
	!Burgers' equation using different numerical methods
	
	use globals
	
	write(*,*) 'Numerical methods for BURGERS EQUATION'
	write(*,*) 'Inviscid (i) or viscid (v) Burgers Equation'
	read(*,*) what_burgers_equation
	
	if (what_burgers_equation .EQ. 'v') then
		what_method=7
	elseif (what_burgers_equation .EQ. 'i') then
		write(*,*) 'What numerical method should be used'
		write(*,*) ' 1 natural finite difference method'
		write(*,*) ' 2 upwind method'
		write(*,*) ' 3 Lax Friedrich method'
		write(*,*) ' 4 Godunovs method'
		write(*,*) ' 5 Lax Wendroff method'
		write(*,*) ' 6 High resolution - slope limiter method'
		read(*,*) what_method
		if (what_method .NE. 1 .AND. what_method .NE. 2 .AND. what_method .NE. 3 .AND. &
			& what_method .NE. 4 .AND. what_method .NE. 5 .AND. what_method .NE. 6) then
			write(*,*) 'bad value for what_method'
			stop
		endif
	else
		write(*,*) 'bad value for what_burgers_equation'
		stop
	endif
	
	call init ()
	
	call init_data ()
	
	call calculation ()
	
	call put_out ()
	
end program burgers_equation

subroutine init ()
	use globals
	implicit none
	
	open(unit=input,file='input.a') !This does not seem to exist
	read(input,*) a
	read(input,*) b
	read(input,*) n
	read(input,*) margin
	read(input,*) courant_number
	read(input,*) t_max
	if (what_method .EQ. 7) read(input,*) D
	close(input)
	
	allocate(x(1:n))
	allocate(h(1:n-1))
	allocate(u(1:n))
	allocate(u_old(1:n))
	allocate(numerical_flux(1:n-1))
	
	n_intv=n-1
	n_intv_i=1.0/n_intv
	
	do i=1,n
		x(i)=(i-1)*(b-a)*n_intv_i
	enddo
	do i=1,n_intv
		h(i)=x(i+1)-x(i)
	enddo
	
	delta_x=(b-a)*n_intv_i
	delta_t=delta_x*courant_number
	
end subroutine init

subroutine init_data ()
	use globals
	implicit none
	write(*,*) 'Initial data'
	write(*,*) '	(h) roof initial data'
	write(*,*) '	(s) shock initial data'
	write(*,*) '	(r) rarefaction initial data'
	write(*,*) '	(f) shock forming initial data'
	write(*,*) '	(c) a complicated example'
	read(*,*) initial_data
	
	if (initial_data .EQ. 'h') then
		l_per=.true.
		quarter=n/4
		half=n/2
		do i=1,quarter
			u(i)=0
		enddo
		do i=n-quarter+1,n
			u(i)=0
		enddo
		do i=quarter+1,half+1
			u(i)=(x(i)-x(quarter+1)) / (x(half+1)-x(quarter+1))
		enddo
		do i=half+1,n-quarter
			u(i)=1-(x(i)-x(half+1))/(x(n-quarter)-x(half+1))
		enddo
	elseif (initial_data .EQ. 's') then
		l_per=.false.
		quarter=n/4
		u_l=1.2
		u_r=0.4
		do i=1,quarter
			u(i)=u_l
		enddo
		do i=quarter+1,n
			u(i)=u_r
		enddo
	elseif (initial_data .EQ. 'r') then
		l_per=.false.
		quarter=n/4
		u_l=1.2
		u_r=0.4
		do i=1,quarter
			u(i)=u_l
		enddo
	elseif (initial_data .EQ. 'f') then
		l_per=.false.
		quarter=n/4
		half=n/2
			do i=1,quarter
			u(i)=0.25
		enddo
		do i=quarter+1,half+1
			u(i)=0.25 - (x(i)-x(quarter+1))*0.25/(x(half+1)-x(quarter+1))
		enddo
		do i=half+2,n
			u(i)=0
		enddo
	elseif (initial_data .EQ. 'c') then
		l_per=.true.
		quarter=n/4
		half=n/2
		do i=1,quarter
			u(i)=0
		enddo
		do i=quarter+1,half+1
			u(i)=1
		enddo
		do i=half+2,n
			u(i)=0
		enddo
	else
		write(*,*) 'bad value for initial_data'
		stop
	endif
end subroutine init_data

subroutine put_out ()
	use globals
	implicit none
	
	open(unit=output,file='output.xls',form='formatted',position='rewind',status='unknown')
	do i=1,n
		write(output,*) u(i)
	enddo
end subroutine put_out

subroutine calculation ()
	use globals
	implicit none
	do time=0,t_max,delta_t !time,t_max,delta_t must be integers but are defined as real - must investigate usage
		if (first_time .EQV. .false.) then
			if (l_per) then
				u(1:margin)=u(n-margin-margin+1:n-margin)
				u(n-margin+1:n)=u(margin+1:margin+margin)
			else
				if (what_method .NE. 6) then
					u(1)=u(2)
					u(n)=u(n-1)
				endif
				if (what_method .EQ. 6) then
					u(1)=u(3)
					u(2)=u(3)
					u(n)=u(n-2)
					u(n-1)=u(n-2)
				endif
			endif
		endif
		first_time=.false.
		
		u_old=u
		
		if (what_method .EQ. 1) call finite_difference ()
		if (what_method .EQ. 2) call upwind ()
		if (what_method .EQ. 3) call lax_friedrich ()
		if (what_method .EQ. 4) call godunov ()
		if (what_method .EQ. 5) call lax_wendroff ()
		if (what_method .EQ. 6) call high_resolution ()
		if (what_method .EQ. 7) call parabolic ()
	enddo
end subroutine calculation

subroutine finite_difference ()
	use globals
	implicit none
	do i=2,n
		u(i)=u_old(i)-courant_number*u_old(i)*(u_old(i)-u_old(i-1))
	enddo
end subroutine finite_difference

subroutine upwind ()
	use globals
	implicit none
	do i=2,n
		u(i)=u_old(i)-courant_number*(0.5*(u_old(i)**2)-0.5*(u_old(i-1)**2))
	enddo
end subroutine upwind

subroutine lax_friedrich
	use globals
	implicit none
	do i=2,n-1
		u(i)=0.5*(u_old(i-1) + u_old(i+1) - 0.5*courant_number*(f(u_old(i+1))-f(u_old(i-1))))
	enddo
end subroutine lax_friedrich

subroutine godunov
	use globals
	implicit none
	
	!calculating numerical flux
	do i=1,n-1
		if (f_der(u_old(i)) .GE. 0 .AND. f_der(u_old(i+1)) .GE. 0) u_star=u_old(i)
		if (f_der(u_old(i)) .LE. 0 .AND. f_der(u_old(i+1)) .LE. 0) u_star=u_old(i+1)
		if (f_der(u_old(i)) .GE. 0 .AND. f_der(u_old(i+1)) .LE. 0 .AND. &
			& (f(u_old(i+1)-f(u_old(i)))/(u_old(i+1)-u_old(i))) .GT. 0) u_star=u_old(i)
		if (f_der(u_old(i)) .GE. 0 .AND. f_der(u_old(i+1)) .LE. 0 .AND. &
			& (f(u_old(i+1)-f(u_old(i)))/(u_old(i+1)-u_old(i))) .LT. 0) u_star=u_old(i+1)	
		if (f_der(u_old(i)) .LT. 0 .AND. f_der(u_old(i+1)) .GT. 0) u_star=0
		numerical_flux(i)=f(u_star)
	enddo
	
	do i=2,n-1
		u(i)=u_old(i)-courant_number*(numerical_flux(i)-numerical_flux(i-1))
	enddo
end subroutine godunov

subroutine parabolic
	use globals
	implicit none
	do i=2,n-1
		f_minus=0.5*(f(u_old(i))+f(u_old(i-1)))
		f_plus=0.5*(f(u_old(i))+f(u_old(i+1)))
		u(i)=u_old(i)+courant_number*(D*(u_old(i+1)-2*u_old(i)+u_old(i-1))/delta_x - &
			& (f_plus-f_minus)) 
!			& (f(u_old(i)-f(u_old(i-1))))
	enddo
end subroutine parabolic

subroutine lax_wendroff ()
	use globals
	implicit none
	do i=2,n-1
		u(i)=u_old(i)-0.5*courant_number*(f(u_old(i+1))-f(u_old(i-1)))+0.5*courant_number**2* &
			& (f_der(0.5*(u_old(i)+u_old(i+1)))*(f(u_old(i+1))-f(u_old(i))) - &
			& f_der(0.5*(u_old(i)+u_old(i-1)))*(f(u_old(i))-f(u_old(i-1))))
	enddo
end subroutine lax_wendroff

subroutine high_resolution ()
	use globals
	implicit none
	
	!calculation numerical flux
	do i=2,n-1
		u1=u_old(i+1)-u_old(i)
		u2=u_old(i)-u_old(i-1)
		if ( abs(u1) .LT. abs(u2) .AND. (u1*u2) .GT. 0) minmod=u1
		if ( abs(u2) .LT. abs(u1) .AND. (u1*u2) .GT. 0) minmod=u2
		if ( (u1*u2) .LE. 0) minmod=0
		sigma=1/delta_x*minmod
		if ( abs(u_old(i+1)-u_old(i)) .GT. 0.001) then
			a_roof=(f(u_old(i+1))-f(u_old(i)))/(u_old(i+1)-u_old(i))
		else
			a_roof=0
		endif
		numerical_flux(i)=f(u_old(i)) + 0.5*a_roof*(1-courant_number*a_roof)*delta_x*sigma
	enddo
	
	do i=3,n-2
		u(i)= u_old(i) -courant_number*(numerical_flux(i)-numerical_flux(i-1))
	enddo
end subroutine high_resolution