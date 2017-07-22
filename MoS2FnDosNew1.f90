!PERIODIC GB as DISLOCATION ARRAYS
!阉割版三带二硫化钼gap计算
!consider only the nearest neighbour Mo atoms near a periodic grain boundary
include "libdf1.f90"
module MatOps

contains

   subroutine insertblock(a, b, matrix, block, x, y) !插入矩阵块
      implicit none
      integer :: a, b, x, y
      complex :: matrix(a*b, a*b), block(a, a)
      integer :: i, j, m, n
!~   matrix((x-1)*a+1:x*a,(y-1)*a+1:y*a)=block
      do i = (x - 1)*a + 1, x*a
         do j = (y - 1)*a + 1, y*a
            m = i - (x - 1)*a
            n = j - (y - 1)*a
            matrix(i, j) = block(m, n)
         end do
      end do
   end subroutine
   subroutine addblock(a, b1, b2, matrix, block, x, y) !插入矩阵块
      !我操,忘了只能插入方阵
      implicit none
      integer :: a, b1, b2, x, y
      complex :: matrix(a*b1, a*b2), block(a, a)
      integer :: i, j, m, n
!~   matrix((x-1)*a+1:x*a,(y-1)*a+1:y*a)=matrix((x-1)*a+1:x*a,(y-1)*a+1:y*a)+block
      do i = (x - 1)*a + 1, x*a
         do j = (y - 1)*a + 1, y*a
            m = i - (x - 1)*a
            n = j - (y - 1)*a
            matrix(i, j) = matrix(i, j) + block(m, n)
         end do
      end do
   end subroutine
!读取格式：
!a,---1stN-------,2ndN,3rdN

end module

module parm
   integer, parameter :: ns = 33, nk = ns !*ns
end module

module modperiod
   use parm
   integer, parameter :: ndos = 255, ndv = 5, nld = 2, nrd = 2
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

   function expr(k)
      implicit none
      real :: k
      real, parameter :: pi1 = -3.1416*2
      complex :: expr
      expr = exp((0, pi1)*k)
      return
   end function

   subroutine surf(g0, n, hm0, v, z)
      implicit none
      integer :: n, i, j
      complex :: g0(n, n), hm1(n, n), v(n, n), vd(n, n), sf(n, n), ivg(n, n), z, hm0(n, n)
      vd = conjg(transpose(v))
      do i = 1, 100
         hm1 = 0
         sf = matmul(v, g0)
         sf = matmul(sf, vd)
         hm1 = hm0 + sf
         g0 = greene(n, hm1, z)
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
      hcpl = hcpl + conjg(transpose(hcpl))
   end subroutine
   !单个方向的周期性:不同lv方向各有一个耦合矩阵,所以以上只包含一个方向（需要逐个计算）
   !注意：周期耦合的矩阵只需要给出a->b 分块的分量而不需要同时包含其对称

   function hamk(hamkp, n, kk) result(hmk)
      implicit none
      integer :: n
      real :: kk
      complex :: hamkp(n, n), hmk(n, n)
      real, parameter :: pi1 = -3.1416*2
      hmk = hamkp*exp((0, pi1)*kk)
      hmk = hmk + conjg(transpose(hmk))
      return
   end function

   function greene(n, ham, z) result(gf)
      !生成一定k下的格林函数
      !k确定时相当于在单胞上加了一个修正项;只需考虑单胞内情况
      integer :: n
      complex :: z, ham(n, n), gf(n, n), ivg(n, n)

      ivg = eye(n, z)
      ivg = ivg - ham
      gf = 0
      gf = cInverse(n, ivg)
      return
   end function

   !确定k->计算lead自能->加入到device蛤密顿H(k)中->1/(z-Hk)积分得到device格林函数->格林函数对k积分：实空间格林函数

   subroutine dosMo(z, dosi, hmdv, hmdpd, hmrd, hmrdpd, hmrdcp, hmrdcp1, hmrdcp2, hmld, hmldpd, hmldcp, hmldcp1, hmldcp2,&
   &hmlcp, hmlcp1, hmlcp2, hmrcp, hmrcp1, hmrcp2, hdfct, kk, kdos)

      implicit none
      real, parameter :: pi1 = -3.1416*2
      integer :: i, j, k, tmps(3), ttmp
      real :: kk(nk, 2), kp, dosi, kdos(nk)
      complex :: sfedv(ndv*3, ndv*3), selfeg(ndv*3, ndv*3)
      complex :: hmdv(ndv*3, ndv*3), hmdpd(ndv*3, ndv*3)
      complex :: hmrd(nrd*3, nrd*3), hmrdpd(nrd*3, nrd*3), hmrdcp(nrd*3, nrd*3)&
      &, hmrdcp1(nrd*3, nrd*3), hmrdcp2(nrd*3, nrd*3)
      complex :: hmld(nld*3, nld*3), hmldpd(nld*3, nld*3), hmldcp(nld*3, nld*3)&
      &, hmldcp1(nld*3, nld*3), hmldcp2(nld*3, nld*3)
      complex :: hmlcp(ndv*3, nld*3), hmlcp1(ndv*3, nld*3), hmlcp2(ndv*3, nld*3)
      complex :: hmrcp(ndv*3, nrd*3), hmrcp1(ndv*3, nrd*3), hmrcp2(ndv*3, nrd*3)
      complex :: hdfct(3*ndv, 3*ndv)

      complex :: hamkr(nrd*3, nrd*3), hamkl(nrd*3, nrd*3), hkpdr(3*nrd, 3*nrd), hkpdl(3*nld, 3*nld)
      complex :: hkpdv(3*ndv, 3*ndv), hmkdv(3*ndv, 3*ndv)
      complex :: hrcp1(3*nrd, 3*ndv), hlcp1(3*nld, 3*ndv)

      complex :: gdv(ndv*3, ndv*3), gdvk(ndv*3, ndv*3), gr(nrd*3, nrd*3), gl(nld*3, nld*3)

      complex :: gredf(ndv*3, ndv*3), ivg(ndv*3, ndv*3), z, eppx
      !需要的矩阵：内部3个（左右D）；周期couple3个（左右D）；连接线4个（LL,LD,RR,RD)
      !所以一共使用10个矩阵……比patch方法多得多了……………………

      !装载k空间采样(根据superlattice的倒易矢量：正交性)
      !kspace采样的范围：设晶格矢量为1，则k分量范围从-pi到pi
!~ print *,kk
      gdv = 0
      write (*, *) z
      do i = 1, nk
         kp = kk(i, 2)
!~           write(*,*) kp
         call hampd(hmrdpd, nrd*3, kp, hkpdr)
         call hampd(hmldpd, nld*3, kp, hkpdl)
         hamkr = hmrd + hkpdr
         hamkl = hmld + hkpdl
         gr = greene(nrd*3, hamkr, z)
         gl = greene(nld*3, hamkl, z)

         hmrdcp = hmrdcp + hmrdcp1*expr(kp) + hmrdcp2*expr(-kp) !+hmrdcp2*expr(-kp)
         hmldcp = hmldcp + hmldcp1*expr(kp) + hmldcp2*expr(-kp)
         call surf(gr, nrd*3, hamkr, hmrdcp, z)
         call surf(gl, nld*3, hamkl, hmldcp, z) !给定k得到lead自能

         !得出左右lead耦合的作用（自能）
         hmlcp = hmlcp + hmlcp1*expr(kp) + hmlcp2*expr(-kp)
         hmrcp = hmrcp + hmrcp1*expr(kp) + hmrcp2*expr(-kp) !expr(-kp)
         !Coupling with lead terms in the next transverse cell

         sfedv = 0
         hlcp1 = conjg(transpose(hmlcp))
         sfedv = sfedv + matmul(hmlcp, matmul(gl, hlcp1))
         hrcp1 = conjg(transpose(hmrcp))
         sfedv = sfedv + matmul(hmrcp, matmul(gr, hrcp1))
         call hampd(hmdpd, ndv*3, kp, hkpdv)
         hmkdv = hkpdv + sfedv + hmdv
         gdvk = greene(ndv*3, hmkdv, z)
         kdos(i) = 0
         do j = 1, ndv*3
            kdos(i) = kdos(i) - imag(gdvk(j, j))
         end do
         gdv = gdv + gdvk/nk
      end do
      write (*, *) "offset", z
      dosi = 0
!~   do i=1,ndv
!~                 if(rmved(i)==0) dosi=dosi-imag(gredf(i,i))
!~         end do

   end subroutine

end module

!自能terms是否具有周期性?

program main
   use modperiod
   use MatOps
   implicit none

   integer :: i, j, k
   integer :: tn1(ndv), ttmp, rmved(ndv)
   real :: dosi(ndos), kdos(ndos, nk)
   real :: kk(nk, 2)
   complex :: z(ndos)
   real, parameter :: t0 = -0.184, t1 = 0.401, t2 = 0.507, t11 = 0.218, t12 = 0.338, t22 = 0.057, s3 = sqrt(3.0)
   complex, parameter :: h01(3, 3) = reshape([t0, t1, t2, -t1, t11, t12, t2, -t12, t22], [3, 3]),&
   &h03(3, 3) = reshape([t0, 0.5*t1 - s3*t2/2.0, -s3*t1/2.0-0.5*t2,&
   &-0.5*t1 - s3*t2/2.0, 0.25*t11+0.75*t22, -s3*0.25*t11-t12+s3*t22/4.0,&
   &s3*t1/2.0-0.5*t2, -s3*t11/4.0+t12+s3*t22/4.0, 0.75*t11+0.25*t22], [3, 3]),&
   &h02(3, 3) = reshape([t0, 0.5*t1 + s3*t2/2.0, s3*t1/2.0-0.5*t2,&
   &-0.5*t1 + s3*t2/2.0, 0.25*t11+0.75*t22, s3*0.25*t11-t12-s3*t22/4.0,&
   &-s3*t1/2.0-0.5*t2, s3*t11/4.0+t12-s3*t22/4.0, 0.75*t11+0.25*t22], [3, 3])
   complex, parameter :: h1(3, 3) = transpose(h01), h2(3, 3) = transpose(h02), h3(3, 3) = transpose(h03)
   complex, parameter :: ep1(3, 3) = reshape([1.046, 0.0, 0.0, 0.0, 2.104, 0.0, 0.0, 0.0, 2.104], [3, 3])
   complex :: hmdv(ndv*3, ndv*3), hmdpd(ndv*3, ndv*3)
   complex :: hmrd(nrd*3, nrd*3), hmrdpd(nrd*3, nrd*3), hmrdcp(nrd*3, nrd*3), hmrdcp1(nrd*3, nrd*3), hmrdcp2(nrd*3, nrd*3)
   complex :: hmld(nld*3, nld*3), hmldpd(nld*3, nld*3), hmldcp(nld*3, nld*3), hmldcp1(nld*3, nld*3), hmldcp2(nld*3, nld*3)
   complex :: hmlcp(ndv*3, nld*3), hmlcp1(ndv*3, nld*3), hmlcp2(ndv*3, nld*3)
   complex :: hmrcp(ndv*3, nrd*3), hmrcp1(ndv*3, nld*3), hmrcp2(ndv*3, nld*3)
   complex :: hdfct(3*ndv, 3*ndv), hblks(-3:3, 3, 3)
   real, parameter :: ehgh = 2.13, elow = -0.56

   open (unit=233, file='fmatrix')

   hblks = 0
   hblks(1, :, :) = h1
   hblks(2, :, :) = h2
   hblks(3, :, :) = h3
   hblks(-1, :, :) = transpose(h1)
   hblks(-2, :, :) = transpose(h2)
   hblks(-3, :, :) = transpose(h3)

!~ print *,h1
!~ stop

   !
   do k = 1, 1
      hmdv = 0
      do i = 1, ndv
         read (233, *) ttmp, (tn1(j), j=1, ndv)
         do j = 1, ndv
            call addblock(3, ndv, ndv, hmdv, hblks(tn1(j), :, :), i, j)
         end do
         call addblock(3, ndv, ndv, hmdv, ep1, i, i)
      end do
      hmdpd = 0
      do i = 1, ndv
         read (233, *) ttmp, (tn1(j), j=1, ndv)
         do j = 1, ndv
            call addblock(3, ndv, ndv, hmdpd, hblks(tn1(j), :, :), i, j)
         end do
      end do

      hmrd = 0
      do i = 1, nrd
         read (233, *) ttmp, (tn1(j), j=1, nrd)
         do j = 1, nrd
            call addblock(3, nrd, nrd, hmrd, hblks(tn1(j), :, :), i, j)
         end do
         call addblock(3, nrd, nrd, hmrd, ep1, i, i)
      end do
      hmrdpd = 0
      do i = 1, nrd
         read (233, *) ttmp, (tn1(j), j=1, nrd)
         do j = 1, nrd
            call addblock(3, nrd, nrd, hmrdpd, hblks(tn1(j), :, :), i, j)
         end do
      end do
      hmrdcp = 0
      do i = 1, nrd
         read (233, *) ttmp, (tn1(j), j=1, nrd)
         do j = 1, nrd
            call addblock(3, nrd, nrd, hmrdcp, hblks(tn1(j), :, :), i, j)
         end do
      end do
      hmrdcp1 = 0
      do i = 1, nrd
         read (233, *) ttmp, (tn1(j), j=1, nrd)
         do j = 1, nrd
            call addblock(3, nrd, nrd, hmrdcp1, hblks(tn1(j), :, :), i, j)
         end do
      end do
      hmrdcp2 = 0
      do i = 1, nrd
         read (233, *) ttmp, (tn1(j), j=1, nrd)
         do j = 1, nrd
            call addblock(3, nrd, nrd, hmrdcp2, hblks(tn1(j), :, :), i, j)
         end do
      end do

      !lOAD THE LEFT LEAD HAMILTONIANS
      hmld = 0
      do i = 1, nld
         read (233, *) ttmp, (tn1(j), j=1, nld)
         do j = 1, nld
            call addblock(3, nld, nld, hmld, hblks(tn1(j), :, :), i, j)
         end do
         call addblock(3, nld, nld, hmld, ep1, i, i)
      end do
      hmldpd = 0
      do i = 1, nld
         read (233, *) ttmp, (tn1(j), j=1, nld)
         do j = 1, nld
            call addblock(3, nld, nld, hmldpd, hblks(tn1(j), :, :), i, j)
         end do
      end do
      hmldcp = 0
      do i = 1, nld
         read (233, *) ttmp, (tn1(j), j=1, nld)
         do j = 1, nld
            call addblock(3, nld, nld, hmldcp, hblks(tn1(j), :, :), i, j)
         end do
      end do
      hmldcp1 = 0
      do i = 1, nld
         read (233, *) ttmp, (tn1(j), j=1, nld)
         do j = 1, nld
            call addblock(3, nld, nld, hmldcp1, hblks(tn1(j), :, :), i, j)
         end do
      end do
      hmldcp2 = 0
      do i = 1, nld
         read (233, *) ttmp, (tn1(j), j=1, nld)
         do j = 1, nld
            call addblock(3, nld, nld, hmldcp2, hblks(tn1(j), :, :), i, j)
         end do
      end do

      hmrcp = 0
      do i = 1, ndv
         read (233, *) ttmp, (tn1(j), j=1, nrd)
         do j = 1, nrd
            call addblock(3, ndv, nrd, hmrcp, hblks(tn1(j), :, :), i, j)
         end do
      end do
      hmrcp1 = 0
      do i = 1, ndv
         read (233, *) ttmp, (tn1(j), j=1, nrd)
         do j = 1, nrd
            call addblock(3, ndv, nrd, hmrcp1, hblks(tn1(j), :, :), i, j)
         end do
      end do
      hmrcp2 = 0
      do i = 1, ndv
         read (233, *) ttmp, (tn1(j), j=1, nrd)
         do j = 1, nrd
            call addblock(3, ndv, nrd, hmrcp2, hblks(tn1(j), :, :), i, j)
         end do
      end do

      hmlcp = 0
      do i = 1, ndv
         read (233, *) ttmp, (tn1(j), j=1, nld)
         do j = 1, nld
            call addblock(3, ndv, nld, hmlcp, hblks(tn1(j), :, :), i, j)
         end do
      end do
      hmlcp1 = 0
      do i = 1, ndv
         read (233, *) ttmp, (tn1(j), j=1, nld)
         do j = 1, nld
            call addblock(3, ndv, nld, hmlcp1, hblks(tn1(j), :, :), i, j)
         end do
      end do
      hmlcp2 = 0
      do i = 1, ndv
         read (233, *) ttmp, (tn1(j), j=1, nld)
         do j = 1, nld
            call addblock(3, ndv, nld, hmlcp2, hblks(tn1(j), :, :), i, j)
         end do
      end do

!~                 end do
   end do

!~         rewind(233)
   open (unit=235, file='kspace.csv')
   do i = 1, nk
      read (235, *) kk(i, 1), kk(i, 2)
   end do
!~    open (unit=234, file='defdev')
!~    hdfct = 0
!~    do i = 1, ndv
!~       read (234, *) ttmp, (tn1(j), j=1, 6), rmved(i)
!~       do j = 1, 6
!~          if (tn1(j) /= 0) call addblock(3, ndv, hdfct, h1, i, tn1(j))
!~       end do
!~    end do
   z = (0, 0)

   do i = 1, ndos
      dosi(i) = 0.0
   end do
!$OMP         PARALLEL DO
   do i = 1, ndos
!~                 write(*,*) i,z
      z(i) = elow + i*(ehgh - elow)/real(ndos) + (0.0, 0.01)
      call dosMo(z(i), dosi(i), hmdv, hmdpd, hmrd, hmrdpd, hmrdcp, hmrdcp1, hmrdcp2, hmld, hmldpd, hmldcp, hmldcp1, hmldcp2&
      &, hmlcp, hmlcp1, hmlcp2, hmrcp, hmrcp1, hmrcp2, hdfct, kk, kdos(i, :))
   end do
   open (unit=49, file='densos.csv')
   open (unit=45, file='densosK.csv')
   open (unit=473, file='dokk.csv')
   do i = 1, ndos
      write (49, *) real(z(i)), ',', dosi(i)
      write (45, *) (kdos(i, j), ',', j=1, nk)
      do j = 1, nk
         write (473, *) real(z(i)), ',', j, ',', kdos(i, j)
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
