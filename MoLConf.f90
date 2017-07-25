module param
	
   integer, parameter :: row =5, crow = 4
   integer, parameter :: ndv = (crow-1)*row+(row+1)/2
   integer,parameter :: layer=1
	integer,parameter :: nld=layer*row,nrd=layer*row
	integer,parameter :: deg=6
	contains
	function indx(ro,cl)
		implicit none
		integer :: i,j,indx,ro,cl
		indx=0
		do i=1,ro-1
			indx=indx+crow-1+mod(i,2)
		end do
		indx=indx+cl
		return
	end function
	function shift(xy,d) result(nxy)
		implicit none
		integer :: xy(2),d(2),dx,dy,nxy(2)
		dx=d(1)
		dy=d(2)
		nxy(1)=xy(1)+dx
		nxy(2)=xy(2)+dy
		return
	end function

	function shiftdev(n,xy,d) result(dxy)
		implicit none
		integer :: n,xy(n,2),d(2),dxy(n,2)
		integer :: i,j,k
		dxy=0
		do i=1,n
			dxy(i,:)=shift(xy(i,:),d)
		end do
		return
	end function

	function hcoupling(dev1,dev2,n1,n2,hcp) result(ham)
		implicit none
		integer :: n1,n2,dist(2)
		integer :: dev1(n1,2),dev2(n2,2),ham(n1,n2)
		integer :: i,j,k,hcp(-100:100,-100:100)
		ham=0
		do i=1,n1
			do j=1,n2
				dist=(dev1(i,:)-dev2(j,:))
				ham(i,j)=hcp(dist(1),dist(2))
!~ 				write(*,*) i,j,dist(1),dist(2),ham(i,j)
			end do
		end do
		return
	end function
!写尼玛什么函数,果断空间换时间
end module
include "networklibstandard.f90"

program main
	use param
	implicit none
	integer :: i,j,k,lns,lna,lnb,va
	integer :: ttmp,defect(ndv,ndv),removed(ndv),reshp(0:ndv)
	integer :: dev(ndv,2),ld(nld,2),rd(nrd,2)
	integer :: devp(ndv,2),ldp(nld,2),rdp(nrd,2)
	integer :: vsh(13,2),dist(2)
	integer :: hdev(0:ndv,0:ndv),nep
	integer :: hcp(-100:100,-100:100)=0
	integer,parameter :: hh(-1:1,-1:1)=reshape([0,-2,3,-1,0,1,-3,2,0],[3,3])
	integer :: hmdv(ndv,ndv)
	character(7) :: inputs
	hcp(-1:1,-1:1)=hh

	do i=0,ndv 
		reshp(i)=i
	end do 
nep=ndv
	
	
	do i=1,row
		do j=1,crow
			dev(indx(i,j),2)=i-1
			dev(indx(i,j),1)=j-1
		end do
	end do
		hmdv=hcoupling(dev,dev,ndv,ndv,hcp)
removed=0
!~ 	print *,indx(3,2),dev(6,:)
!~ 	stop
	201 read(*,*) inputs
	select case(inputs)
		case('rmv')
			read(*,*) lns
			do i=1,ndv
				defect(i,lns)=0
				defect(lns,i)=0
			end do
			removed(lns)=1
			reshp(lns)=0
			do i=lns+1,ndv
				reshp(i)=reshp(i)-1
			end do
			nep=nep-1
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

	

	202 open(unit=233,file='thedev')
	hdev=0
	
	do i=1,ndv
		do j=1,ndv
			hdev(reshp(i),reshp(j))=hmdv(i,j)
		end do
	end do 
	do i=1,nep
		write(233,*) i, (hdev(i,j),j=1,ndv)
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
