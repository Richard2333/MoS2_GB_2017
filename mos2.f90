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
  matrix((x-1)*a+1:x*a,(y-1)*a+1:y*a)=matrix((x-1)*a+1:x*a,(y-1)*a+1:y*a)+block
!~   do i=(x-1)*a+1,x*a
!~     do j=(y-1)*a+1,y*a
!~       m=i-(x-1)*a
!~ 			n=j-(y-1)*a
!~ 			matrix(i,j)=block(m,n)
!~ 		end do
!~ 	end do
end subroutine
!读取格式：
!a,---1stN-------,2ndN,3rdN



end module
