!PERIODIC GB as DISLOCATION ARRAYS
!阉割版三带二硫化钼gap计算
!consider only the nearest neighbour Mo atoms near a periodic grain boundary
include "libdf1.f90"
module MatOps

contains

subroutine insertblock(a,b,matrix,block,x,y)!插入矩阵块
	implicit none
	integer :: a,b,x,y
	complex :: matrix(a*b,a*b),block(a,a)
  integer :: i,j,m,n
  matrix((x-1)*a+1:x*a,(y-1)*a+1:y*a)=block
!~   do i=(x-1)*a+1,x*a
!~     do j=(y-1)*a+1,y*a
!~       m=i-(x-1)*a
!~ 			n=j-(y-1)*a
!~ 			matrix(i,j)=block(m,n)
!~ 		end do
!~ 	end do
end subroutine
subroutine addblock(a,b,matrix,block,x,y)!插入矩阵块
	implicit none
	integer :: a,b,x,y
	complex :: matrix(a*b,a*b),block(a,a)
  integer :: i,j,m,n
!~   matrix((x-1)*a+1:x*a,(y-1)*a+1:y*a)=matrix((x-1)*a+1:x*a,(y-1)*a+1:y*a)+block
  do i=(x-1)*a+1,x*a
    do j=(y-1)*a+1,y*a
      m=i-(x-1)*a
			n=j-(y-1)*a
			matrix(i,j)=matrix(i,j)+block(m,n)
		end do
	end do
end subroutine
!读取格式：
!a,---1stN-------,2ndN,3rdN



end module

module parm
  integer, parameter :: nslice = 3, row = 4, crow = 8, ns = 33, nk = ns!*ns
end module
module modq
  use parm
contains
  function eye(n, z) result(dig)
    integer :: n, i
    complex :: z, dig(n, n)
    dig = 0
    do i = 1, n
      dig(i, i) = z
    end do
    !  return
  end function
  function f(k) !d1=(0,1)
    implicit none
    complex :: f
    real :: k(2)
    real, parameter :: pi1 = 3.1415926*2
    f = 0
    f = f + exp((0, pi1)*k(2))
    f = f + exp((0, pi1)*(0.866*k(1) - 0.5*k(2)))
    f = f + exp((0, pi1)*(-0.866*k(1) - 0.5*k(2)))
    f = -2.9*f
    return
  end function f
  function cInverse(n, a) result(ra)
    integer::n, lda, ipiv(n), info, lwork
    complex ::a(n, n), ra(n, n), work(n)
    ra = a
    lwork = n
    lda = n
    call cgetrf(n, n, ra, lda, ipiv, info)
    if (info /= 0) write (0, *) 'Error occured in cgetrf!', info
    call cgetri(n, ra, lda, ipiv, work, lwork, info)
    if (info /= 0) write (0, *) 'Error occured in cgetri!', info
  end function
  subroutine kspace1(k)
    implicit none
    integer, parameter :: nk = 6*4**nslice
    real :: k(nk, 2)
    real :: dlx, dly, my, mx
    integer :: i, j, nx, ny, idx

    !nslice=2
    ny = 2**nslice
    mx = 2/(sqrt(3.0)*3.0)
    dlx = mx/ny

    nx = ny
    idx = 0
    my = -mx*sqrt(3.0)/2
    dly = -my/ny
    !write(*,*) dlx/dly
    do i = 1, ny
      do j = -nx, nx
        idx = idx + 1
        k(idx, 1) = j*dlx/2
        k(idx, 2) = my + dly/3 + mod(j + nx + 1, 2)*dly/3
      end do
      nx = nx + 1
      my = my + dly
    end do
    do i = 1, idx
      k(i + idx, 1) = k(i, 1)
      k(i + idx, 2) = -k(i, 2)
    end do
    idx = idx*2

    write (*, *) idx
  end subroutine kspace1

  subroutine kspace(k) !Generate the k-space samples with monkhorst scheme
    implicit none

    real :: k(nk, 2), k1(2), k2(2)
    real :: dlx, dly, my, mx, nx, ny
    integer :: i, j, idx, idy, ind

    !nslice=2

    k1(2) = 1.0/3
    k1(1) = sqrt(3.0)/3
    k2(1) = k1(1)
    k2(2) = -k1(2)
    k = 0
    ind = 0
    do i = 0, ns - 1
      do j = 1, ns
        idx = (2*(i + 1) - ns - 1)
        idy = (2*j - ns - 1)
        write (*, *) idx, idy, idx + idy
        nx = 1.0*idx/(2*ns)
        ny = 1.0*idy/(2*ns)
        k(i*ns + j, :) = nx*k1 + ny*k2
        !if(k(i*ns+j,2)>k2(2)*0.5) k(i*ns+j,2)=k(i*ns+j,2)-k2(2)
        !if(k(i*ns+j,2)<-k2(2)*0.5) k(i*ns+j,2)=k(i*ns+j,2)+k2(2)
      end do
    end do
    write (*, *) 'innn', ind

  end subroutine kspace

  function greenz(k, nk, xy, z, n, r1, r2)
    implicit none
    real :: k(nk, 2), xy(n, 3), dx, dy, kk(2)
    integer :: nk, i, n, r1, r2
    complex :: z, greenz, four1, ab
    real, parameter :: pi1 = 3.1415926*2
    !complex,external :: f
    dx = xy(r1, 1) - xy(r2, 1)
    dy = xy(r1, 2) - xy(r2, 2)
    greenz = 0
    do i = 1, nk
      kk(1) = k(i, 1)
      kk(2) = k(i, 2)

      if (xy(r1, 3) == xy(r2, 3)) ab = z
      if (xy(r1, 3) == 1 .and. xy(r2, 3) == 0) ab = f(kk)
      if (xy(r2, 3) == 1 .and. xy(r1, 3) == 0) ab = conjg(f(kk))
      four1 = exp((0, pi1)*(kk(1)*dx + kk(2)*dy))
      !write(*,*) four1
      greenz = greenz + four1*ab/(z**2 - abs(f(kk))**2)
      !write(55,*) i,greenz,(k(i,1)*dx+k(i,2)*dy),k(i,1),k(i,2)
    end do
    !write(*,*) 'a',greenz
    greenz = greenz/nk
    !write(*,*) 'b',greenz

    return
  end function greenz

  subroutine genHa(latce, ham, n)
    !On-site potential,Vi=omg*(1-sumLDOS)=omg*localdense
    !use mod1
    implicit none
    integer :: latce(n, 3), n, i, j
    complex :: ham(n, n)
    real, parameter :: tS = 2.9
    ham = 0
    do i = 1, n
      do j = 1, 3
        if (latce(i, j) /= 0) then
          ham(i, latce(i, j)) = (-1, 0)*tS
        end if
      end do

    end do

  end subroutine genHa

  !对一个superlattice根据周期性生成动量k的哈密顿量
  !k的红外截断给出superlattice的周期性(n个superlattice为一个周期),紫外截断给出superlattice的大小

  function ivmato(mat, n0, n1, n)
    implicit none
    integer :: n0, n1, n
    complex :: mat(n, n), ra(n1 - n0 + 1, n1 - n0 + 1), rb(n1 - n0 + 1, n1 - n0 + 1), ivmato(n, n)
    !complex,external cInverse
    ra = 0
    ra = mat(n0:n1, n0:n1)
    rb = cInverse(n1 - n0 + 1, ra)
    ivmato = 0
    ivmato(n0:n1, n0:n1) = rb

    return
  end function ivmato

  function degree(lattice, n, shell)
    implicit none
    integer :: i, j, n, lattice(n, 3), degree(n), shell(n)
    degree = 0
    do i = 1, n
      do j = 1, 3
        if (lattice(i, j) /= 0 .and. shell(lattice(i, j)) == 0) degree(i) = degree(i) + 1
      end do
    end do
    return
  end function degree

  subroutine reshap_mat(mat, n, resp, mat1)
    implicit none
    integer :: n, resp(n), i, j
    complex :: mat(n, n), mat1(n, n)
    do i = 1, n
      do j = 1, n
        mat1(resp(i), resp(j)) = mat(i, j)
      end do
    end do
  end subroutine reshap_mat

  subroutine selfeng(latce, greens, z, n, selfe)
    implicit none
    integer :: n
    integer :: latce(n, 3), i, j
    !integer :: respmat(n,n)
    complex :: ham(n, n), z, greens(n, n)
    complex :: ivgret(n, n), ivg(n, n), selfe(n, n), gen(n, n), sfe
    call genHa(latce, ham, n)
    ivgret = 0
    do i = 1, n
      ivgret(i, i) = z
    end do

    ivgret = ivgret - ham
    gen = greens
    ivg = cInverse(n, greens)
    selfe = ivgret - ivg


    do i = 1, n
      do j = 1, n
        write (37, *) i, ',', j, ',', real(gen(i, j)), ',', imag(gen(i, j))
      end do

    end do
  end subroutine selfeng

  subroutine reshape_lattice(lattice, resp, n, latc)
    implicit none
    integer :: n
    integer :: lattice(n, 3), i, j, shell(n), degs(n), resp(n), cnt, latc(n, 3), icl, n0
    !integer :: respmat(n,n)
    !~         integer,parameter :: nk=6*4**nslice

    shell = 0
    cnt = 0
    icl = 0
    n0 = n

    degs = degree(lattice, n, shell)
    do i = 1, n
      if (degs(i) /= 3 .and. shell(i) == 0) then
        cnt = cnt + 1
        resp(i) = cnt
        shell(i) = 1
        icl = icl + 1
      end if
    end do
    !n0=n0-icl

    do i = 1, n
      if (shell(i) == 0) then
        cnt = cnt + 1
        resp(i) = cnt
        shell(i) = 1
      end if
    end do

    latc = 0
    do i = 1, n
      do j = 1, 3
        if (lattice(i, j) /= 0) then
          latc(resp(i), j) = resp(lattice(i, j))
          !对应：i->resp(i)
          !lattice(i,j)->resp(lattice(i,j)
          !于是有：i,lattice(i,j)=1 -> resp(i),resp(lattice(i,j))=1
        end if
      end do
    end do
  end subroutine reshape_lattice

  subroutine gensfe(lattice, k, xy, z, selfe)
    integer, parameter :: n = row*crow
    real :: k(nk, 2), xy(n, 3)
    integer :: lattice(n, 3), resp(n), i, j, latc(n, 3)
    complex :: selfe(n, n), z, greens(n, n), green1(n, n)
    call kspace(k)
    do i = 1, n
      do j = 1, n
        greens(i, j) = greenz(k, nk, xy, z, n, i, j)
      end do
    end do
    !~         call reshap_mat(greens,n,resp,green1)
    call selfeng(lattice, greens, z, n, selfe)
  end subroutine gensfe

end module modq

module modperiod
  use modq
  integer ,parameter :: ndos=255,ndv=18,nld=4,nrd=4
contains

	function expr(k)
		implicit none
		real :: k
		real,parameter :: pi1=-3.1416*2
		complex :: expr
		expr=exp((0,pi1)*k)
		return
	end function

  subroutine surf(g0, n, hm0, v, z)
    implicit none
    integer :: n, i, j
    complex :: g0(n, n), hm1(n, n), v(n, n), vd(n, n), sf(n, n), ivg(n, n), z,hm0(n,n)
    vd = conjg(transpose(v))
    do i=1,32
    	hm1=0
      sf = matmul(v, g0)
      sf = matmul(sf, vd)
      hm1 = hm0 + sf
      g0=greene(n,hm1,z)
    end do

  end subroutine
  !能否用k得出一条线的GF之后再进行recursive:可以!因为recursive的自能只包括两点格林函数s
  subroutine hampd(hamcp, n, kk, hcpl)
    implicit none
    integer :: n, i, j
    complex :: hamcp(n, n), hcpl(n, n)
    real :: kk
    real, parameter :: pi1 = -3.1416*2
    hcpl = hamcp*exp((0, pi1)*kk)
    hcpl=hcpl+conjg(transpose(hcpl))
  end subroutine
  !单个方向的周期性:不同lv方向各有一个耦合矩阵,所以以上只包含一个方向（需要逐个计算）
  !注意：周期耦合的矩阵只需要给出a->b 分块的分量而不需要同时包含其对称

	function hamk(hamkp,n,kk) result(hmk)
		implicit none
		integer :: n
		real :: kk
		complex :: hamkp(n,n),hmk(n,n)
		real,parameter :: pi1=-3.1416*2
		hmk=hamkp*exp((0,pi1)*kk)
		hmk=hmk+conjg(transpose(hmk))
		return
	end function

  function greene(n,ham,z) result(gf)
    !生成一定k下的格林函数
    !k确定时相当于在单胞上加了一个修正项;只需考虑单胞内情况
    integer :: n
    complex :: z,ham(n,n),gf(n,n),ivg(n,n)

    ivg=eye(n,z)
    ivg=ivg-ham
    gf=0
    gf=cInverse(n,ivg)
    return
  end function


  !确定k->计算lead自能->加入到device蛤密顿H(k)中->1/(z-Hk)积分得到device格林函数->格林函数对k积分：实空间格林函数

subroutine dosMo(z,dosi,hmdv,hmdpd,hmrd,hmrdpd,hmrdcp,hmrdcp1,hmrdcp2,hmld,hmldpd,hmldcp,hmldcp1,hmldcp2,&
&hmlcp,hmlcp1,hmlcp2,hmrcp,hmrcp1,hmrcp2,hdfct,kk,rmved,kdos)
  
  implicit none
  real,parameter :: pi1=-3.1416*2
  integer :: i,j,k,tmps(3),ttmp
  real :: kk(nk,2),kp,dosi,kdos(nk)
  complex :: sfedv(ndv*3,ndv*3),selfeg(ndv*3,ndv*3)
	complex :: hmdv(ndv*3,ndv*3),hmdpd(ndv*3,ndv*3)
	complex :: hmrd(nrd*3,nrd*3),hmrdpd(nrd*3,nrd*3),hmrdcp(nrd*3,nrd*3)&
	&,hmrdcp1(nrd*3,nrd*3),hmrdcp2(nrd*3,nrd*3)
	complex :: hmld(nld*3,nld*3),hmldpd(nld*3,nld*3),hmldcp(nld*3,nld*3)&
	&,hmldcp1(nld*3,nld*3),hmldcp2(nld*3,nld*3)
	complex :: hmlcp(ndv*3,nld*3),hmlcp1(ndv*3,nld*3),hmlcp2(ndv*3,nld*3)
	complex :: hmrcp(ndv*3,nrd*3),hmrcp1(ndv*3,nrd*3),hmrcp2(ndv*3,nrd*3)
	complex :: hdfct(3*ndv,3*ndv)

	complex :: hamkr(nrd*3,nrd*3),hamkl(nrd*3,nrd*3),hkpdr(3*nrd,3*nrd),hkpdl(3*nld,3*nld)
	complex :: hkpdv(3*ndv,3*ndv),hmkdv(3*ndv,3*ndv)
	complex :: hrcp1(3*nrd,3*ndv),hlcp1(3*nld,3*ndv)
	
	complex :: gdv(ndv*3,ndv*3),gdvk(ndv*3,ndv*3),gr(nrd*3,nrd*3),gl(nld*3,nld*3)
	integer :: rmved(ndv)
  complex,parameter :: t=(-2.8,0)
	complex :: gredf(ndv*3,ndv*3),ivg(ndv*3,ndv*3),z,eppx
  !需要的矩阵：内部3个（左右D）；周期couple3个（左右D）；连接线4个（LL,LD,RR,RD)
  !所以一共使用10个矩阵……比patch方法多得多了……………………

  !装载k空间采样(根据superlattice的倒易矢量：正交性)
  !kspace采样的范围：设晶格矢量为1，则k分量范围从-pi到pi
!~ print *,kk
  gdv=0
  write(*,*) z
  do i=1,nk
  	kp=kk(i,2)
!~   	write(*,*) kp
    call hampd(hmrdpd,nrd*3,kp,hkpdr)
    call hampd(hmldpd,nld*3,kp,hkpdl)
    hamkr=hmrd+hkpdr
    hamkl=hmld+hkpdl
    gr=greene(nrd*3,hamkr,z)
    gl=greene(nld*3,hamkl,z)

	eppx=exp((0,pi1)*kp)
	
		hmrdcp=hmrdcp+hmrdcp1*eppx+hmrdcp2*conjg(eppx)!+hmrdcp2*conjg(eppx)
		hmldcp=hmldcp+hmldcp1*eppx+hmldcp2*conjg(eppx)
    call surf(gr,nrd*3,hamkr,hmrdcp,z)
    call surf(gl,nld*3,hamkl,hmldcp,z)!给定k得到lead自能

		!得出左右lead耦合的作用（自能）
		hmlcp=hmlcp+hmlcp1*eppx+hmlcp2*conjg(eppx)
		hmrcp=hmrcp+hmrcp1*eppx+hmrcp2*conjg(eppx)!conjg(eppx)
		!Coupling with lead terms in the next transverse cell
    
    sfedv=0
    hlcp1=conjg(transpose(hmlcp))
    sfedv=sfedv+matmul(hmlcp,matmul(gl,hlcp1))
    hrcp1=conjg(transpose(hmrcp))
    sfedv=sfedv+matmul(hmrcp,matmul(gr,hrcp1))
    call hampd(hmdpd,ndv*3,kp,hkpdv)
    hmkdv=hkpdv+sfedv+hmdv
		gdvk=greene(ndv*3,hmkdv,z)
		kdos(i)=0
		do j=1,ndv*3
			kdos(i)=kdos(i)-imag(gdvk(j,j))
		end do
    gdv=gdv+gdvk/nk
  end do
    write(*,*) "offset",z

!~   ivg=cInverse(ndv*3,gdv)
!~ 	selfeg=eye(ndv*3,z)
!~ 	selfeg=selfeg-ivg-hmdv


!~   selfeg=cInverse(ndv,gdv)

!~   stop


!~   ivg=eye(ndv*3,z)
!~   ivg=ivg-hdfct-selfeg
!~   gredf=0
!~   gredf=cInverse(ndv*3,ivg)
!~   !生成格点格林函数（缺陷核心区）
  dosi=0
!~   do i=1,ndv
!~ 		if(rmved(i)==0) dosi=dosi-imag(gredf(i,i))
!~ 	end do

end subroutine



end module

!自能terms是否具有周期性?



program main
	use modperiod
	use MatOps
	implicit none

	integer :: i,j,k
	complex,parameter :: t=(-2.8,0)
	integer :: tn1(ndv),ttmp,rmved(ndv)
	real :: dosi(ndos),kdos(ndos,nk)
	real :: kk(nk,2)
	complex :: z(ndos)
	real,parameter :: t0=-0.184,t1=0.401,t2=0.507,t11=0.218,t12=0.338,t22=0.057,s3=sqrt(3.0)
  complex,parameter :: h1(3,3)=reshape([t0,t1,t2,-t1,t11,t12,t2,-t12,t22],[3,3]),&
  &h3(3,3)=reshape([t0,0.5*t1-s3*t2/2.0,-s3*t1/2.0-0.5*t2,&
  &-0.5*t1-s3*t2/2.0,0.25*t11+0.75*t22,-s3*0.25*t11-t12+s3*t22/4.0,&
  &s3*t1/2.0-0.5*t2,-s3*t11/4.0+t12+s3*t22/4.0,0.75*t11+0.25*t22],[3,3]),&
  &h2(3,3)=reshape([t0,0.5*t1+s3*t2/2.0,s3*t1/2.0-0.5*t2,&
  &-0.5*t1+s3*t2/2.0,0.25*t11+0.75*t22,s3*0.25*t11-t12-s3*t22/4.0,&
  &-s3*t1/2.0-0.5*t2,s3*t11/4.0+t12-s3*t22/4.0,0.75*t11+0.25*t22],[3,3])
  complex,parameter :: ep1(3,3)=reshape([1.046,0.0,0.0,0.0,2.104,0.0,0.0,0.0,2.104],[3,3])
	complex :: hmdv(ndv*3,ndv*3),hmdpd(ndv*3,ndv*3)
	complex :: hmrd(nrd*3,nrd*3),hmrdpd(nrd*3,nrd*3),hmrdcp(nrd*3,nrd*3),hmrdcp1(nrd*3,nrd*3),hmrdcp2(nrd*3,nrd*3)
	complex :: hmld(nld*3,nld*3),hmldpd(nld*3,nld*3),hmldcp(nld*3,nld*3),hmldcp1(nld*3,nld*3),hmldcp2(nld*3,nld*3)
	complex :: hmlcp(ndv*3,nld*3),hmlcp1(ndv*3,nld*3),hmlcp2(ndv*3,nld*3)
	complex :: hmrcp(ndv*3,nrd*3),hmrcp1(ndv*3,nld*3),hmrcp2(ndv*3,nld*3)
	complex :: hdfct(3*ndv,3*ndv),hblks(-3:3,3,3)
	real,parameter :: ehgh=2.13,elow=-0.56
	
	open(unit=233,file='fmatrix')

	hblks=0
	hblks(1,:,:)=h1
	hblks(2,:,:)=h2
	hblks(3,:,:)=h3
	hblks(-1,:,:)=transpose(h1)
	hblks(-2,:,:)=transpose(h2)
	hblks(-3,:,:)=transpose(h3)
	
!~ print *,h1
!~ stop

	!
!~ do k=1,1
		hmdv=0	
		do i=1,ndv
			read(233,*) ttmp,(tn1(j),j=1,ndv)
			do j=1,ndv
				call addblock(3,ndv,hmdv,hblks(tn1(j),:,:),i,j)
			end do
			call addblock(3,ndv,hmdv,ep1,i,i)
		end do
!~ 		do i=1,ndv
!~ 			write(77,*) (real(hmdv(i,j)),j=1,ndv)
!~ 		end do
!~ 		stop
		hmdpd=0
		do i=1,ndv
			read(233,*) ttmp,(tn1(j),j=1,ndv)
			do j=1,ndv
				call addblock(3,ndv,hmdpd,hblks(tn1(j),:,:),i,j)
			end do
		end do

		hmrd=0
		do i=1,nrd
			read(233,*) ttmp,(tn1(j),j=1,nrd)
			do j=1,nrd
				call addblock(3,nrd,hmrd,hblks(tn1(j),:,:),i,j)
			end do
			call addblock(3,nrd,hmrd,ep1,i,i)
		end do 
		hmrdpd=0
		do i=1,nrd
			read(233,*) ttmp,(tn1(j),j=1,nrd)
			do j=1,nrd
				call addblock(3,nrd,hmrdpd,hblks(tn1(j),:,:),i,j)
			end do
		end do 
		hmrdcp=0
		do i=1,nrd
			read(233,*) ttmp,(tn1(j),j=1,nrd)
			do j=1,nrd
				call addblock(3,nrd,hmrdcp,hblks(tn1(j),:,:),i,j)
			end do
		end do
		hmrdcp1=0
		do i=1,nrd
			read(233,*) ttmp,(tn1(j),j=1,nrd)
			do j=1,nrd
				call addblock(3,nrd,hmrdcp1,hblks(tn1(j),:,:),i,j)
			end do
		end do
		hmrdcp2=0
		do i=1,nrd
			read(233,*) ttmp,(tn1(j),j=1,nrd)
			do j=1,nrd
				call addblock(3,nrd,hmrdcp2,hblks(tn1(j),:,:),i,j)
			end do
		end do

		!lOAD THE LEFT LEAD HAMILTONIANS
		hmld=0
		do i=1,nld
			read(233,*) ttmp,(tn1(j),j=1,nld)
			do j=1,nld
				call addblock(3,nld,hmld,hblks(tn1(j),:,:),i,j)
			end do
			call addblock(3,nld,hmld,ep1,i,i)
		end do 
		hmldpd=0
		do i=1,nld
			read(233,*) ttmp,(tn1(j),j=1,nld)
			do j=1,nld
				call addblock(3,nld,hmldpd,hblks(tn1(j),:,:),i,j)
			end do
		end do 
		hmldcp=0
		do i=1,nld
			read(233,*) ttmp,(tn1(j),j=1,nld)
			do j=1,nld
				call addblock(3,nld,hmldcp,hblks(tn1(j),:,:),i,j)
			end do
		end do
		hmldcp1=0
		do i=1,nld
			read(233,*) ttmp,(tn1(j),j=1,nld)
			do j=1,nld
				call addblock(3,nld,hmldcp1,hblks(tn1(j),:,:),i,j)
			end do
		end do
		hmldcp2=0
		do i=1,nld
			read(233,*) ttmp,(tn1(j),j=1,nld)
			do j=1,nld
				call addblock(3,nld,hmldcp2,hblks(tn1(j),:,:),i,j)
			end do
		end do

		hmrcp=0
		do i=1,ndv
			read(233,*) ttmp,(tn1(j),j=1,nrd)
			do j=1,nrd
				call addblock(3,nld,hmrcp,hblks(tn1(j),:,:),i,j)
			end do
		end do
		hmrcp1=0
		do i=1,ndv
			read(233,*) ttmp,(tn1(j),j=1,nrd)
			do j=1,nrd
				call addblock(3,nld,hmrcp1,hblks(tn1(j),:,:),i,j)
			end do
		end do
		hmrcp2=0
		do i=1,ndv
			read(233,*) ttmp,(tn1(j),j=1,nrd)
			do j=1,nrd
				call addblock(3,nld,hmrcp2,hblks(tn1(j),:,:),i,j)
			end do
		end do
		
		hmlcp=0
		do i=1,ndv
			read(233,*) ttmp,(tn1(j),j=1,nld)
			do j=1,nld
				call addblock(3,nld,hmlcp,hblks(tn1(j),:,:),i,j)
			end do
		end do
		hmlcp1=0
		do i=1,ndv
			read(233,*) ttmp,(tn1(j),j=1,nld)
			do j=1,nld
				call addblock(3,nld,hmlcp1,hblks(tn1(j),:,:),i,j)
			end do
		end do
		hmlcp2=0
		do i=1,ndv
			read(233,*) ttmp,(tn1(j),j=1,nld)
			do j=1,nld
				call addblock(3,nld,hmlcp2,hblks(tn1(j),:,:),i,j)
			end do
		end do
		
!~ 		end do
!~ end do 

	
!~ 	rewind(233)
  open(unit=235,file='kspace.csv')
  do i=1,nk
    read(235,*) kk(i,1),kk(i,2)
  end do
!~   rewind(235)
	open(unit=234,file='defdev')
  hdfct=0
  do i=1,ndv
    read(234,*) ttmp,(tn1(j),j=1,6),rmved(i)
    do j=1,6
				if(tn1(j)/=0) call addblock(3,ndv,hdfct,h1,i,tn1(j))
			end do
  end do
!~   rewind(234)



	z=(0,0)

	do i=1,ndos
		dosi(i)=0.0
	end do
	!$OMP PARALLEL DO
	do i=1,ndos
!~ 		write(*,*) i,z
	z(i)=elow+i*(ehgh-elow)/real(ndos)+(0.0,0.01)
		call  dosMo(z(i),dosi(i),hmdv,hmdpd,hmrd,hmrdpd,hmrdcp,hmrdcp1,hmrdcp2,hmld,hmldpd,hmldcp,hmldcp1,hmldcp2&
		&,hmlcp,hmlcp1,hmlcp2,hmrcp,hmrcp1,hmrcp2,hdfct,kk,rmved,kdos(i,:))
	end do
	open(unit=49,file='densos.csv')
	open(unit=45,file='densosK.csv')
	open(unit=473,file='dokk.csv')
	do i=1,ndos
		write(49,*) real(z(i)),',',dosi(i)
		write(45,*) (kdos(i,j),',',j=1,nk)
		do j=1,nk
			write(473,*) real(z(i)),',',j,',',kdos(i,j)
		end do
	end do


end program main

subroutine coordi(i, crow, row, x, y)
  implicit none
  integer :: i, crow, row, nx, ny
  real :: x, y, ab
  nx = i/crow
  ny = mod(i, crow)
  if (ny == 0) then
    nx = nx - 1
    ny = crow
  end if

  x = 1.732*ny/2
  y = 1.5*(nx) + (-1)**(nx + ny + 1)*0.25
  if (mod(nx, 2) == mod(ny, 2)) ab = 0
  if (mod(nx, 2) /= mod(ny, 2)) ab = 1
end subroutine coordi
