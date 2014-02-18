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
			write(*.*) 'bad value for what_method')
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
	
	open(unit=input,file='input.a')
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
		h(i)-x(i+1)-x(i)
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
	do i=n-quarter_1,n
		u(i)=0
	enddo
	do i=quarter+1,half+1
		u(i)=(x(i)-x(quarter+1))/(x(half+1)-x(quarter+1))
	enddo
	do ihalf+1,n-quarter
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
	do i-quarter+1,half+1
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
	do i=half2,n
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
	do time-0,t_max,delta_t
		if (first_time .EQ. .false.) then
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
		fist_time=.false.
		
		u_old=u
		
		if (what_method .EQ. 1) call finite_difference ()
		if (what_method .EQ. 2) call upwind ()
		if (what_method .EQ. 3) call lax_friedrich ()
		if (what_method .EQ. 4) call godunov ()
		if (what_method .EQ. 5) call lax_wendroff ()
		if (what_method .EQ. 6) call high_resolution ()
		if (what_method .EQ. 7) call parabolic ()
	enddo
end subroutine calculation	!this is on page 48 of the paper