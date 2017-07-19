!PERIODIC GB as DISLOCATION ARRAYS
include "libdf1.f90"
module parm
  integer, parameter :: nslice = 3, row = 4, crow = 8, ns = 234, nk = ns!*ns
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
contains
  subroutine surf(g0, n, hm0, v, z)
    implicit none
    integer :: n, i, j
    complex :: g0(n, n), hm1(n, n), v(n, n), vd(n, n), sf(n, n), ivg(n, n), z,hm0(n,n)
    vd = conjg(transpose(v))
    do i=1,400
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



end module

!自能terms是否具有周期性?

program main
  use modperiod
  implicit none
  integer :: i,j,k,tmps(3),ttmp
  integer,parameter :: ndv=22,nrd=6,nld=6,n=ndv
  real :: kk(nk,2),kp
  complex :: sfedv(ndv,ndv),selfeg(ndv,ndv)
  complex :: hmd(ndv,ndv),hmdpd(ndv,ndv)
  complex :: hmrd(nrd,nrd),hmrdpd(nrd,nrd),hmrdcp(nrd,nrd)
  complex :: hmld(nld,nld),hmldpd(nld,nld),hmldcp(nld,nld)
	complex :: hmlcp(ndv,nld),hmrcp(ndv,nrd),hrcp1(nrd,ndv),hlcp1(nld,ndv)
	complex :: hmkdv(ndv,ndv),hkpdv(ndv,ndv),hkpdr(nrd,nrd),hkpdl(nld,nld),hamkr(nrd,nrd),hamkl(nld,nld)
	complex :: gdv(ndv,ndv),gdvk(ndv,ndv),gr(nrd,nrd),gl(nld,nld)
	real :: dfct(ndv,ndv)
  complex,parameter :: t=(-2.8,0),z=(0,0.01)
	complex :: gredf(ndv,ndv),hdfct(ndv,ndv),ivg(ndv,ndv)
	integer :: rmved(ndv)
  !需要的矩阵：内部3个（左右D）；周期couple3个（左右D）；连接线4个（LL,LD,RR,RD)
  !所以一共使用10个矩阵……比patch方法多得多了……………………
  open(unit=233,file='hmatrix')
    !device
      hmd=0
      do i=1,ndv
          read(233,*) ttmp,(tmps(j),j=1,3)
          do j=1,3
            if(tmps(j)/=0) hmd(i,tmps(j))=t
          end do
      end do
      hmdpd=0
      do i=1,ndv
        read(233,*) ttmp,(tmps(j),j=1,3)
        do j=1,3
          if(tmps(j)/=0) hmdpd(i,tmps(j))=t
        end do
      end do
    !right lead
      hmrd=0
      do i=1,nrd
        read(233,*) ttmp,(tmps(j),j=1,3)
        do j=1,3
          if(tmps(j)/=0) hmrd(i,tmps(j))=t
        end do
      end do
      hmrdpd=0
      do i=1,nrd
        read(233,*) ttmp,(tmps(j),j=1,3)
        do j=1,3
          if(tmps(j)/=0) hmrdpd(i,tmps(j))=t
        end do
      end do
      hmrdcp=0
      do i=1,nrd
        read(233,*) ttmp,(tmps(j),j=1,3)
        do j=1,3
          if(tmps(j)/=0) hmrdcp(i,tmps(j))=t
        end do
      end do
    !left lead
      hmld=0
      do i=1,nld
        read(233,*) ttmp,(tmps(j),j=1,3)
        do j=1,3
          if(tmps(j)/=0) hmld(i,tmps(j))=t
        end do
      end do
      hmldpd=0
      do i=1,nld
        read(233,*) ttmp,(tmps(j),j=1,3)
        do j=1,3
          if(tmps(j)/=0) hmldpd(i,tmps(j))=t
        end do
      end do
      hmldcp=0
      do i=1,nld
        read(233,*) ttmp,(tmps(j),j=1,3)
        do j=1,3
          if(tmps(j)/=0) hmldcp(i,tmps(j))=t
        end do
      end do
    !couple
      hmlcp=0
      do i=1,nld
        read(233,*) ttmp,(tmps(j),j=1,3)
        do j=1,3
          if(tmps(j)/=0) hmlcp(tmps(j),i)=t
        end do
      end do
      hmrcp=0
      do i=1,nrd
        read(233,*) ttmp,(tmps(j),j=1,3)
        do j=1,3
          if(tmps(j)/=0) hmrcp(tmps(j),i)=t
        end do
      end do
  !载入晶格矩阵

  open(unit=235,file='kspace.csv')
  do i=1,nk
    read(235,*) kk(i,1),kk(i,2)
  end do
  !装载k空间采样(根据superlattice的倒易矢量：正交性)
  !kspace采样的范围：设晶格矢量为1，则k分量范围从-pi到pi
!~ print *,kk
  gdv=0
  do i=1,nk
  	kp=kk(i,2)
  	write(*,*) kp
    call hampd(hmrdpd,nrd,kp,hkpdr)
    call hampd(hmldpd,nld,kp,hkpdl)
    hamkr=hmrd+hkpdr
    hamkl=hmld+hkpdl
    gr=greene(nrd,hamkr,z)
    gl=greene(nld,hamkl,z)

    call surf(gr,nrd,hamkr,hmrdcp,z)
    call surf(gl,nld,hamkl,hmldcp,z)

    !给定k得到lead自能
    sfedv=0
    hlcp1=conjg(transpose(hmlcp))
    sfedv=sfedv+matmul(hmlcp,matmul(gl,hlcp1))
    hrcp1=conjg(transpose(hmrcp))
    sfedv=sfedv+matmul(hmrcp,matmul(gr,hrcp1))
    call hampd(hmdpd,ndv,kp,hkpdv)
    hmkdv=hkpdv+sfedv+hmd
		gdvk=greene(ndv,hmkdv,z)
    gdv=gdv+gdvk/nk
  end do
  ivg=cInverse(ndv,gdv)
	selfeg=eye(ndv,z)
	selfeg=selfeg-ivg-hmd
  

!~   selfeg=cInverse(ndv,gdv)


  open(unit=234,file='defdev')
  hdfct=0
  do i=1,ndv
    read(234,*) ttmp,(tmps(j),j=1,3),rmved(i)
    do j=1,3
      if(tmps(j)/=0) hdfct(i,tmps(j))=t
    end do
  end do
  ivg=eye(ndv,z)
  ivg=ivg-hdfct-selfeg
  gredf=cInverse(ndv,ivg)

	open(unit=44,file='selfenergy.0710')
	open(unit=45,file='devicegreen.0710')
	 do i=1,ndv
	 	if(rmved(i)==0) then
  	do j=1,ndv
  		if(rmved(j)==0)then
				write(44,*) i,',',j,',',real(selfeg(i,j)),',',imag(selfeg(i,j))
				write(45,*) i,',',j,',',real(gredf(i,j)),',',imag(gredf(i,j))
			end  if
  	end do
  	write(47,*) i,',',real(gredf(i,i)),',',imag(gredf(i,i))
  	end if
  end do
  
  !生成格点格林函数（缺陷核心区）


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
