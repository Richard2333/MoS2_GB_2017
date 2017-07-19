module param
	
   integer, parameter :: row =7, crow = 9
   integer, parameter :: sca = row*crow
	integer,parameter :: deg=3
end module
include "networklibstandard.f90"
program main
	use param
	use biolibs_fixedscale
	implicit none
	integer :: thenet(sca,deg),ttmp
	integer :: i,j,k,ab,la,lb,lns
	integer :: resp(sca),nremv,removed(sca)
	
	real :: x0,y0,x1,y1,nl,ny,zz,nx,m
	character(5) :: inpt
!~ 	open(unit=19,file='latpos')
	do i=1,sca
		resp(i)=i
	end do
	thenet=0
   do i = 0, row - 1
      do j = 2, crow
         call pedge(thenet, sca, deg, i*crow + j, i*crow + j - 1) !:行内连接
      end do
      !call pedge(thenet,sca,deg,i*crow+1,i*crow+crow)
   end do
	
   do i = 1, row - 1
      do j = 1, crow
         if ((i)*crow + j <= sca) then
            if (mod(i, 2) == 1 .and. mod(j, 2) == 1) then
               call pedge(thenet, sca, deg, i*crow + j, (i - 1)*crow + j)
            end if
            if (mod(i, 2) == 0 .and. mod(j, 2) == 0) then
               call pedge(thenet, sca, deg, i*crow + j, (i - 1)*crow + j)
            end if
         end if
      end do
		end do




   !call deledge(thenet,sca,deg,sca,sca-1)
   !all deledge(thenet,sca,deg,1,2)


	!call disorderzbound(thenet,row,crow,2,0.1)

   open (unit=11, file='latticeb')
   open (unit=13, file='latview1.csv')
	!write the crystal lattice
	nremv=0
	removed=0
	201 read (*,*) inpt
	select case(inpt)
		case('del')
			read(*,*) la,lb
			call deledge(thenet,sca,deg,la,lb)
			go to 201
		case('rmv')
			read(*,*) lns
			call remove(thenet,sca,lns)
			nremv=nremv+1
			do i=lns+1,sca
				resp(i)=resp(i)-1
			end do
			removed(lns)=1
			go to 201
		case('lnk')
			read(*,*) la,lb
			call pedge(thenet,sca,deg,la,lb)
			go to 201
		case default
			go to 202
	end select
	
	
	202 do i=1,sca
!~ 		do j=1,deg
			if(removed(i)/=1) write(11,*) resp(i),(resp(thenet(i,j)),j=1,deg)!(thenet(i,j),',',j=1,deg)
!~ 		end do
	end do
	do i=1,sca
		call coordi(i,crow,row,x0,y0)
!~ 		write(18,*) x0,',',y0,',',i
		do j=1,deg
		if(thenet(i,j)/=0) then
			call coordi(thenet(i,j),crow,row,x1,y1)
			write(13,*) x0,',',y0,',',x1,',',y1
		end if
		
		end do
		
	end do

	do i=1,sca
		m=0
		do j=1,deg
			if(thenet(i,j)/=0) m=1
		end do
		if(m/=0) call coordi(i,crow,row,x0,y0)

 		if(removed(i)/=1) write(19,*) x0,y0!,mod(i+1,2)
	end do
	
end program main

subroutine disorderzbound(thenet,row,crow,zl,ds)
	use biolibs_fixedscale
	implicit none
	integer :: thenet(row*crow,3)
	integer :: row,crow,zl,i,j,k,sca
	real :: ds,s
	call sr1and()
	do i=0,zl-1
		do j=1,crow
			call radm(s)
			if(s<=ds) then
				call remove(thenet,sca,i*crow+j)
				write(*,*) i*crow+j
			end if
		end do 
	end do 

	do i=row-zl,row-1
		do j=1,crow
			call radm(s)
			if(s<=ds) then
				call remove(thenet,sca,i*crow+j)
				write(*,*) i*crow+j
			end if
		end do 
	end do 
end subroutine disorderzbound





subroutine remove(thenet,sca,n)
	use biolibs_fixedscale
	implicit none
	integer :: thenet(sca,3)
	integer :: sca,n1,n2,i,j,maxdeg,n
	maxdeg=3
	n1=n
	do i=1,3
		n2=thenet(n,i)
		if(n2/=0) then

			call deledge(thenet,sca,maxdeg,n1,n2)
		end if 
	end do 
end subroutine remove

subroutine coordi(i,crow,row,x,y)
	implicit none
	integer :: i,crow,row,nx,ny
	real :: x,y
	nx=i/crow
	ny=mod(i,crow)
	if(ny==0) then
		nx=nx-1
		ny=crow
	end if
	!ny=ny+nx
	x=1.732*ny/2
   y = 1.5*(nx) + (-1)**(nx + ny + 1)*0.25

end subroutine coordi
