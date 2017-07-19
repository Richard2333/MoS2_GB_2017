include "networklibstandard.f90"
module defected


integer,parameter :: deg=3

contains
	
	subroutine def_st_wa(lattice,x,y,a,xy)
		!$acc routine seq

	use biolibs_fixedscale

!x行，每行y个
	implicit none
	integer :: x,y,lattice(x*y,3),i,j,k,a,m,xx,yy
	real :: xy(x*y,3)
	


	xx=mod(a,y)
	if(xx==0) xx=y
	yy=a/y+1
	if(mod(yy,2)==mod(xx,2)) then
		m=a+y
	else 
		m=a-y
	end if
	call deledge(lattice,x*y,deg,a,a-1)
	call deledge(lattice,x*y,deg,m,m+1)
	call pedge(lattice,x*y,deg,a,m+1)
	
	call pedge(lattice,x*y,deg,m,a-1)
	xy(m,2)=(xy(m,2)+xy(a,2))/2
	xy(a,2)=xy(m,2)
	xy(a,1)=xy(a-1,1)+0.2
	xy(m,1)=xy(m+1,1)-0.2
	write(*,*) (xy(a,1)+xy(m,1))/2,xy(a,2)
	end subroutine

	subroutine def_st_wa1(lattice,x,y,a,xy)
		!$acc routine seq

	use biolibs_fixedscale

!x行，每行y个
	implicit none
	integer :: x,y,lattice(x*y,3),i,j,k,a,m,xx,yy
	real :: xy(x*y,3)
	


	xx=mod(a,y)
	if(xx==0) xx=y
	yy=a/y+1
	if(mod(yy,2)==mod(xx,2)) then
		m=a+y
	else 
		m=a-y
	end if
	call deledge(lattice,x*y,deg,a,a+1)
	call deledge(lattice,x*y,deg,m,m-1)
	call pedge(lattice,x*y,deg,a,m-1)
	
	call pedge(lattice,x*y,deg,m,a+1)
	xy(m,2)=(xy(m,2)+xy(a,2))/2
	xy(a,2)=xy(m,2)
	xy(a,1)=xy(a-1,1)+0.2
	xy(m,1)=xy(m+1,1)-0.2
	write(*,*) (xy(a,1)+xy(m,1))/2,xy(a,2)
	end subroutine
	

	
	subroutine rmv_atom(lattice,x,y,a)
	use biolibs_fixedscale
		implicit none
		integer :: lattice(x*y,3),a,x,y,i
		do i=1,3
			call deledge(lattice,x*y,deg,a,lattice(a,i))
		end do
	end subroutine
	
end module
