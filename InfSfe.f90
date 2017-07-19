module InftySfeGreen
integer,parameter :: nslice=6,row=8,crow=7,ns=23
	integer,parameter :: n=row*crow,ndf=25,klm=160
	integer,parameter :: nk=ns**2!nk=6*4**nslice
	complex,parameter :: covg=(0,0.002)
contains
function f(k)!d1=(0,1)
	!$acc routine 
	implicit none
	complex :: f
	real :: k(2)
	real,parameter :: pi1=3.1415926*2
	f=0
	f=f+exp((0,pi1)*k(2))
	f=f+exp((0,pi1)*(0.866*k(1)-0.5*k(2)))
	f=f+exp((0,pi1)*(-0.866*k(1)-0.5*k(2)))
	f=-2.9*f
	return
end function f

subroutine kspace(k)!Generate the k-space samples with monkhorst scheme
	implicit none
	real :: k(nk,2),k1(2),k2(2)
	real :: dlx,dly,my,mx,nx,ny
	integer :: i,j,idx,idy,ind
	k1(2)=1.0/3
	k1(1)=sqrt(3.0)/3
	k2(1)=k1(1)
	k2(2)=-k1(2)
	k=0
	ind=0
	do i=0,ns-1
		do j=1,ns
			idx=(2*(i+1)-ns-1)
			idy=(2*j-ns-1)
			!write(*,*) idx,idy,idx+idy
			nx=1.0*idx/(2*ns)
			ny=1.0*idy/(2*ns)
			k(i*ns+j,:)=nx*k1+ny*k2
			!if(k(i*ns+j,2)>k2(2)*0.5) k(i*ns+j,2)=k(i*ns+j,2)-k2(2)
			!if(k(i*ns+j,2)<-k2(2)*0.5) k(i*ns+j,2)=k(i*ns+j,2)+k2(2)
		end do
	end do
	write(*,*) 'innn',ind
	
end subroutine kspace

subroutine kspace1(k)
		implicit none

		real :: k(nk,2)
		real :: dlx,dly,my,mx
		integer :: i,j,nx,ny,idx

	!nslice=2
	ny=2**nslice
	mx=2/(sqrt(3.0)*3.0)
	dlx=mx/ny
	
	nx=ny
	idx=0
	my=-mx*sqrt(3.0)/2
	dly=-my/ny
	!write(*,*) dlx/dly
	do i=1,ny
		do j=-nx,nx
			idx=idx+1
			k(idx,1)=j*dlx/2
			k(idx,2)=my+dly/3+mod(j+nx+1,2)*dly/3
		end do
		nx=nx+1
		my=my+dly
	end do
	do i=1,idx
		k(i+idx,1)=k(i,1)
		k(i+idx,2)=-k(i,2)
	end do
	idx=idx*2

	write(*,*) idx
end subroutine kspace1



function greenz(k,nk,xy,z,n,r1,r2)
!$acc routine 
	implicit none
	real :: k(nk,2),xy(n,3),dx,dy,kk(2)
	integer :: nk,i,n,r1,r2
	complex :: z,greenz,four1,ab
	real,parameter :: pi1=3.1415926*2
	!complex,external :: f
	dx=xy(r1,1)-xy(r2,1)
	dy=xy(r1,2)-xy(r2,2)
	greenz=0
	do i=1,nk
		kk(1)=k(i,1)
		kk(2)=k(i,2)
		if(xy(r1,3)==xy(r2,3)) ab=z
		if(xy(r1,3)==1 .and. xy(r2,3)==0) ab=f(kk)
		if(xy(r2,3)==1 .and. xy(r1,3)==0) ab=conjg(f(kk))
		four1=exp((0,pi1)*(kk(1)*dx+kk(2)*dy))
		!write(*,*) four1
		greenz=greenz+four1*ab/(z**2-abs(f(kk))**2)
		!write(55,*) i,greenz,(k(i,1)*dx+k(i,2)*dy),k(i,1),k(i,2)
	end do
	!write(*,*) 'a',greenz
	greenz=greenz/nk
	!write(*,*) 'b',greenz
	return
end function greenz

subroutine genHa(latce,ham,n)
	!On-site potential,Vi=omg*(1-sumLDOS)=omg*localdense
	!use mod1
	!$acc routine 
	implicit none
	integer :: latce(n,3),n,i,j
	complex :: ham(n,n)
	real,parameter :: tS=2.9
	ham=0
	do i=1,n
		do j=1,3
			if(latce(i,j)/=0) then
				ham(i,latce(i,j))=(-1,0)*tS
			end if
		end do

	end do

end subroutine genHa
!generate TB hamiltonian with on-site potential
function cInverse(n, a)  result(ra)
!$acc routine 
	integer::n,lda,ipiv(n),info,lwork
	complex ::a(n,n),ra(n,n),work(n)
	ra=a
	lwork=n
	lda=n
	call cgetrf(n, n, ra, lda, ipiv, info)
	if(info/=0) write(0,*) 'Error occured in cgetrf!',info
	call cgetri(n, ra, lda, ipiv, work, lwork, info)
	if(info/=0) write(0,*) 'Error occured in cgetri!',info
	endfunction
	

!generate TB hamiltonian with on-site potential

function ivmato(mat,n0,n1,n)
!$acc routine 
	implicit none
	integer :: n0,n1,n
	complex :: mat(n,n),ra(n1-n0+1,n1-n0+1),rb(n1-n0+1,n1-n0+1),ivmato(n,n)
	!complex,external cInverse
	ra=0
	ra=mat(n0:n1,n0:n1)
	rb=cInverse(n1-n0+1,ra)
	ivmato=0
	ivmato(n0:n1,n0:n1)=rb
!	do i=n0,n1
!		do j=n0,n1
!			ivmato(i,j)=rb(i-n0+1,j-n0+1)
!		end do
!	end do

return
end function ivmato

function degree(lattice,n,shell)
!$acc routine 
	implicit none
	integer :: i,j,n,lattice(n,3),degree(n),shell(n)
	degree=0
	do i=1,n
		do j=1,3
			if(lattice(i,j)/=0 .and. shell(lattice(i,j))==0) degree(i)=degree(i)+1
		end do
	end do
	return
end function degree

subroutine reshap_mat(mat,n,resp,mat1)
!$acc routine 
	implicit none
	integer :: n,resp(n),i,j
	complex :: mat(n,n),mat1(n,n)
	do i=1,n
		do j=1,n
			mat1(resp(i),resp(j))=mat(i,j)
		end do 
	end do
end subroutine reshap_mat

subroutine selfeng(latce,greens,z,n,selfe)
!$acc routine 
	implicit none
	integer :: n
	integer :: latce(n,3),i,j
	!integer :: respmat(n,n)
	!integer,parameter :: nk=6*4**nslice
	complex :: ham(n,n),z,greens(n,n)
	complex :: ivgret(n,n),ivg(n,n),selfe(n,n),gen(n,n),sfe
	call genHa(latce,ham,n)
	ivgret=0
	do i=1,n
		ivgret(i,i)=z
	end do
	
	ivgret=ivgret-ham
	gen=greens
	ivg=cInverse(n,greens)
	selfe=ivgret-ivg
end subroutine selfeng

subroutine reshape_lattice(lattice,resp,n,latc)
!$acc routine 
	implicit none
	integer :: n
	integer :: lattice(n,3),i,j,shell(n),degs(n),resp(n),cnt,latc(n,3),icl,n0
	!integer :: respmat(n,n)
	!integer,parameter :: nk=6*4**nslice
	shell=0
	cnt=0
	icl=0
	n0=n
	degs=degree(lattice,n,shell)
	do i=1,n
		if(degs(i)/=3 .and. shell(i)==0) then
			cnt=cnt+1
			resp(i)=cnt
			shell(i)=1
			icl=icl+1
		end if
	end do
	!n0=n0-icl
	do i=1,n
		if(shell(i)==0) then
			cnt=cnt+1
			resp(i)=cnt
			shell(i)=1
		end if
	end do
	
	latc=0
	do i=1,n
		do j=1,3
			if(lattice(i,j)/=0) then
				latc(resp(i),j)=resp(lattice(i,j))
				!对应：i->resp(i)
				!lattice(i,j)->resp(lattice(i,j)
				!于是有：i,lattice(i,j)=1 -> resp(i),resp(lattice(i,j))=1
			end if
		end do
	end do
end subroutine reshape_lattice


subroutine gensfe(lattice,k,xy,z,selfe,kzctr)
!$acc routine 
	integer,parameter :: n=row*crow!,nk=6*4**nslice
	real :: k(nk,2),xy(n,3)
	integer :: lattice(n,3),resp(n),i,j,latc(n,3),kzctr
	complex :: selfe(n,n),z,greens(n,n),green1(n,n)
	call reshape_lattice(lattice,resp,n,latc)
	
	if(kzctr==1) call kspace1(k)
	if(kzctr==0) call kspace(k)
	
	do i=1,n
		do j=1,n
			greens(i,j)=greenz(k,nk,xy,z,n,i,j)
		end do
	end do
	call reshap_mat(greens,n,resp,green1)
	call selfeng(latc,green1,z,n,selfe)
end subroutine gensfe

function dos(latti1,k,xy,z)
!$acc routine 
	use defected
	implicit none
	integer :: lattice(n,3),i,j,latce(n,3),resp(n),nenergy,latti1(n,3)
	real :: k(nk,2),xy(n,3),xy1(n,3),x1,y1,x2,y2
	complex :: z,ham(n,n),ivg(n,n),green(n,n),ham1(n,n),zz(2*klm+1)
	complex :: selfe(n,n),dos
		lattice=latti1
		selfe=0
		call gensfe(lattice,k,xy,z,selfe,1)
		!call def_st_wa(lattice,row,crow,ndf,xy)
		write (*,*) lattice(49,:)
		!call rmv_atom(lattice,row,crow,ndf)
		ham=0
		call reshape_lattice(lattice,resp,n,latce)
		call genHa(latce,ham,n)
		ivg=0
		do i=1,n
			ivg(i,i)=z
		end do
		!selfe=0
		ivg=ivg-ham-selfe
		green=0
		green=cInverse(n,ivg)
		dos=0
!~ 		!call kspace(k)
!~ 		write(*,*) nenergy
!~ 		do i=1,n
!~ 			do j=1,n
!~ 				green(i,j)=greenz(k,nk,xy,z,n,i,j)
!~ 			end do
!~ 		end do
		do i=1,n
			dos=dos+green(i,i)
		end do
		return
end function dos

end module InftySfeGreen
