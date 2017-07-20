module param
	
   integer, parameter :: row =4, crow = 5
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
				dist=dev2(j,:)-dev1(i,:)
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
	integer :: hdev(ndv,ndv)
	integer :: hcp(-100:100,-100:100)=0
	integer,parameter :: hh(-1:1,-1:1)=reshape([0,-2,3,-1,0,1,-3,2,0],[3,3])
	integer :: hmdv(ndv,ndv),hmdpd(ndv,ndv)
	integer :: hmrd(nrd,nrd),hmrdpd(nrd,nrd),hmrdcp(nrd,nrd),hmrdcp1(nrd,nrd),hmrdcp2(nrd,nrd)
	integer :: hmld(nld,nld),hmldpd(nld,nld),hmldcp(nld,nld),hmldcp1(nld,nld),hmldcp2(nld,nld)
	integer :: hmrcp(ndv,nrd),hmrcp1(ndv,nrd),hmrcp2(ndv,nrd)
	integer :: hmlcp(ndv,nld),hmlcp1(ndv,nld),hmlcp2(ndv,nld)
	hcp(-1:1,-1:1)=hh


	
	vsh(1,:)=[-2,4]
	
	vsh(2,:)=[-2,4]
	vsh(3,:)=[1,0]
	vsh(4,:)=[-1,4]
	vsh(5,:)=[3,-4]

	vsh(6,:)=[-2,4]
	vsh(7,:)=[-1,0]
	vsh(8,:)=[-3,4]
	vsh(9,:)=[1,-4]

	vsh(10,:)=[-2,4]
	vsh(11,:)=[2,-4]

	vsh(12,:)=[-2,4]
	vsh(13,:)=[2,-4]
	
	do i=1,row
		do j=1,crow-1+mod(i,2)
			dev(indx(i,j),2)=i-1
			dev(indx(i,j),1)=j-1-(i-1)/2
		end do
	end do 
!~ 	print *,indx(3,2),dev(6,:)
!~ 	stop
	do i=1,row
		rd(i,:)=shift(dev(indx(i,crow-1+mod(i,2)),:),[1,0])
		ld(i,:)=shift(dev(indx(i,1),:),[-1,0])
!~ 		write(*,*) ld(i,:)
	end do
	!periodic of the device region
	hmdv=hcoupling(dev,dev,ndv,ndv,hcp)
	devp=shiftdev(ndv,dev,vsh(1,:))
	hmdpd=hcoupling(dev,devp,ndv,ndv,hcp)

	
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

	hmrcp=hcoupling(dev,rd,ndv,nrd,hcp)
	rdp=shiftdev(nrd,rd,vsh(10,:))
	hmrcp1=hcoupling(dev,rdp,ndv,nrd,hcp)
	rdp=shiftdev(nrd,rd,vsh(11,:))
	hmrcp2=hcoupling(dev,rdp,ndv,nrd,hcp)

	hmlcp=hcoupling(dev,ld,ndv,nld,hcp)
	ldp=shiftdev(nld,ld,vsh(12,:))
	hmlcp1=hcoupling(dev,ldp,ndv,nld,hcp)
	ldp=shiftdev(nld,ld,vsh(13,:))
	hmlcp2=hcoupling(dev,ldp,ndv,nld,hcp)

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
