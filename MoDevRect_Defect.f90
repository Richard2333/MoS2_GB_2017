module param
	
   integer, parameter :: row =5, crow = 5
   integer, parameter :: ndv = 18
   integer,parameter :: layer=1
	integer,parameter :: nld=5,nrd=4
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
	integer :: i,j,k
	integer :: dev(ndv,2),ld(nld,2),rd(nrd,2)
	integer :: devp(ndv,2),ldp(nld,2),rdp(nrd,2)
	integer :: vsh(13,2),dist(2)
	integer :: hdev(ndv,ndv),ttmp
	integer :: hcp(-100:100,-100:100)=0
	integer,parameter :: hh(-1:1,-1:1)=reshape([0,-2,3,-1,0,1,-3,2,0],[3,3])
	integer :: hmdv(ndv,ndv),hmdpd(ndv,ndv)
	integer :: hmrd(nrd,nrd),hmrdpd(nrd,nrd),hmrdcp(nrd,nrd),hmrdcp1(nrd,nrd),hmrdcp2(nrd,nrd)
	integer :: hmld(nld,nld),hmldpd(nld,nld),hmldcp(nld,nld),hmldcp1(nld,nld),hmldcp2(nld,nld)
	integer :: hmrcp(ndv,nrd),hmrcp1(ndv,nrd),hmrcp2(ndv,nrd)
	integer :: hmlcp(ndv,nld),hmlcp1(ndv,nld),hmlcp2(ndv,nld)
	hcp(-1:1,-1:1)=hh

	

	

	
	open(unit=77,file='inleads.csv')
	do i=1,nrd
		read(77,*) ttmp ,rd(i,1),rd(i,2)
	end do
	read(77,*) vsh(2,1),vsh(2,2)
	read(77,*) vsh(3,1),vsh(3,2)
	do i=1,nld
		read(77,*) ttmp, ld(i,1),ld(i,2)
	end do
	read(77,*) vsh(6,1),vsh(6,2)
	read(77,*) vsh(7,1),vsh(7,2)
	

	
	vsh(4,:)=vsh(2,:)+vsh(3,:)
	vsh(5,:)=vsh(3,:)-vsh(2,:)

	vsh(8,:)=vsh(7,:)+vsh(6,:)
	vsh(9,:)=vsh(7,:)-vsh(6,:)

	vsh(10,:)=vsh(2,:)
	vsh(11,:)=-vsh(2,:)

	vsh(12,:)=vsh(6,:)
	vsh(13,:)=-vsh(6,:)
!~ 	print *,indx(3,2),dev(6,:)
!~ 	stop

	!periodic of the device region
	
	hmdv=0
	hmdpd=0
	
	hmrd=hcoupling(rd,rd,nrd,nrd,hcp)
	rdp=shiftdev(nrd,rd,vsh(2,:))
	hmrdpd=hcoupling(rd,rdp,nrd,nrd,hcp)
	rdp=shiftdev(nrd,rd,vsh(3,:))
	hmrdcp=hcoupling(rd,rdp,nrd,nrd,hcp)
	rdp=shiftdev(nrd,rd,vsh(4,:))
	hmrdcp1=hcoupling(rd,rdp,nrd,nrd,hcp)
	rdp=shiftdev(nrd,rd,vsh(5,:))
	hmrdcp2=hcoupling(rd,rdp,nrd,nrd,hcp)

	hmld=hcoupling(ld,ld,nld,nld,hcp)
	ldp=shiftdev(nld,ld,vsh(6,:))
	hmldpd=hcoupling(ld,ldp,nld,nld,hcp)
	ldp=shiftdev(nld,ld,vsh(7,:))
	hmldcp=hcoupling(ld,ldp,nld,nld,hcp)
	ldp=shiftdev(nld,ld,vsh(8,:))
	hmldcp1=hcoupling(ld,ldp,nld,nld,hcp)
	ldp=shiftdev(nld,ld,vsh(9,:))
	hmldcp2=hcoupling(ld,ldp,nld,nld,hcp)

	hmrcp=0
	hmrcp1=0
	hmrcp2=0

	hmlcp=0
	hmlcp1=0
	hmlcp2=0

	open(unit=233,file='fmatrix')
	
	do i=1,ndv
		write(233,*) i, (hmdv(i,j),j=1,ndv)
	end do
	do i=1,ndv
		write(233,*) i, (hmdpd(i,j),j=1,ndv)
	end do
	do i=1,nrd
		write(233,*) i, (hmrd(i,j),j=1,nrd)
	end do
	do i=1,nrd
		write(233,*) i, (hmrdpd(i,j),j=1,nrd)
	end do
	do i=1,nrd
		write(233,*) i, (hmrdcp(i,j),j=1,nrd)
	end do
	do i=1,nrd
		write(233,*) i, (hmrdcp1(i,j),j=1,nrd)
	end do
	do i=1,nrd
		write(233,*) i, (hmrdcp2(i,j),j=1,nrd)
	end do

	do i=1,nld
		write(233,*) i, (hmld(i,j),j=1,nld)
	end do
	do i=1,nld
		write(233,*) i, (hmldpd(i,j),j=1,nld)
	end do
	do i=1,nld
		write(233,*) i, (hmldcp(i,j),j=1,nld)
	end do
	do i=1,nld
		write(233,*) i, (hmldcp1(i,j),j=1,nld)
	end do
	do i=1,nld
		write(233,*) i, (hmldcp2(i,j),j=1,nld)
	end do

	do i=1,ndv
		write(233,*) i, (hmrcp(i,j),j=1,nrd)
	end do
	do i=1,ndv
		write(233,*) i, (hmrcp1(i,j),j=1,nrd)
	end do
	do i=1,ndv
		write(233,*) i, (hmrcp2(i,j),j=1,nrd)
	end do

	
	do i=1,ndv
		write(233,*) i, (hmlcp(i,j),j=1,nld)
	end do
	do i=1,ndv
		write(233,*) i, (hmlcp1(i,j),j=1,nld)
	end do
	do i=1,ndv
		write(233,*) i, (hmlcp2(i,j),j=1,nld)
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
