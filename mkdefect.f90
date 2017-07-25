module param
	
!~ 	integer,parameter :: row=7,crow=10
	integer,parameter :: sca=20
	integer,parameter :: deg=6
end module
include "networklibstandard.f90"
program main
	use param
	implicit none
	integer :: i,j,k,lns,lna,lnb,va
	integer :: ttmp,defect(sca,sca),removed(sca)
	character(7) :: inputs
	
	open(unit=233,file='fmatrix')
	removed=0 
	do i=1,sca
		read(233,*) ttmp,(defect(i,j),j=1,sca)
	end do 
	201 read(*,*) inputs
	select case(inputs)
		case('rmv')
			read(*,*) lns
			do i=1,sca
				defect(i,lns)=0
				defect(lns,i)=0
			end do
			removed(lns)=1 
			go to 201
		case('lnk')
			read(*,*) lna,lnb,va
			
				defect(lna,lnb)=-va
				defect(lnb,lna)=va
			
			removed(lns)=1 
			go to 201
		case('del')
			read(*,*) lna,lnb
			
				defect(lna,lnb)=0
				defect(lnb,lna)=0
			
			removed(lns)=1 
			go to 201
		case default
			go to 202
	end select

	202 open(unit=234,file='defdev')
	do i=1,sca
		write(234,*) i,(defect(i,j),j=1,sca),removed(i)
	end do 
		 
end program main
