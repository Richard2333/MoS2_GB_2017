
module biolibs_fixedscale
  
contains

  subroutine deledge(net2,nend,maxdeg,n1,n2)!delete edge between two nodes
  	!$acc routine seq

    implicit none
    integer ::nend
    integer  :: maxdeg
    integer ::net2(nend,maxdeg)
    integer :: n1,n2
    integer ::i,j,k
    !allocate(net(nend,maxdeg))

    do i=1,maxdeg
       if(net2(n1,i)==n2) then
          net2(n1,i)=0
          !exit
          write(*,*) 'delx1',n1
       end if
    end do
    !write(*,*) 'zz'
    do i=1,maxdeg
       if(net2(n2,i)==n1) then
          net2(n2,i)=0
          !exit
          write(*,*) 'delx2',n2
       end if
    end do
  end subroutine deledge

  subroutine pedge(net2,sca,maxdeg,n1,n2)
	!$acc routine seq

    implicit none
		integer :: sca
    integer  :: maxdeg
    integer ::	net2(sca,maxdeg)
    integer :: n1,n2,p
    integer ::i,j,k
    !allocate(net(nend,maxdeg))
    p=0
    do i=1,maxdeg
       if(net2(n1,i)==n2) p=1
       if(net2(n1,i)==0 .and. p==0) then
          net2(n1,i)=n2
          p=1

       end if
    end do
    p=0
    do i=1,maxdeg
       if(net2(n2,i)==n1) p=1
       if(net2(n2,i)==0 .and. p==0) then
          net2(n2,i)=n1
          p=1

       end if
    end do
  end subroutine pedge

  subroutine testedg(net,sca,n1,n2,maxdeg,test)

    implicit none
    integer :: sca
    integer :: n1,n2,maxdeg,i,test
    integer :: net(sca,maxdeg)
    test=0
    do i=1,maxdeg
       if (net(n1,i) == n2) then
          test=1
          exit
       end if
    end do
  end subroutine testedg

  subroutine grouptes(net,sca,maxdeg,perco,grop) 	!MAGICAL,DONNOT TOUCH!
    implicit none				
    !integer,parameter :: sca=10000	
    integer  :: maxdeg,sca				
    integer  :: net(sca,maxdeg),perco(sca),grop(sca)
    integer :: 	n1(sca)		
    integer :: 	visited(sca),vis(sca),nextstep(sca)
    integer :: 	i,j,k,p,nex,n,sp
    !s1=net(n1,:)

    n1=0

    visited=0
    vis=0
    nextstep=0

    do sp=1,sca
       if(visited(i)==0 .and. perco(i)/=0) then
          visited(i)=1
          nex=1
          nextstep(1)=i
          grop(i)=n
          do while(nex/=0)
             p=nex
             nex=0
             do j=1,p
                vis(j)=nextstep(j)
                nextstep(j)=0
             end do
             do j=1,p
                do k=1,maxdeg
                   if(perco(net(vis(j),k))/=0 .and. net(vis(j),k)/=0 .and. visited(net(vis(j),k))==0) then
                      nex=nex+1
                      nextstep(nex)=net(vis(j),k)
                      visited(net(vis(j),k))=1
                      grop(net(vis(j),k))=n
                   end if
                end do
             end do
          end do
          n=n+1
       end if
    end do

  end subroutine grouptes

  subroutine gropser(grop,sca,gseries) !第n个group计数
    implicit none
    integer :: grop(sca),gseries(sca)
    integer :: i,j,k,n,sca
    gseries=0
    do i=1,sca
       gseries(grop(i))=gseries(grop(i))+1
    end do
  end subroutine gropser





end module biolibs_fixedscale

  subroutine genperco(perco,sca,p)
    implicit none
    integer :: perco(sca),sca
    real :: p,pr
    integer :: i,j,k
    call sr1and()
    call radm(pr)
    do i=1,sca
       call radm(pr)
       if(pr<p) then
          perco(i)=1
       else 
          perco(i)=0
       end if
    end do
  end subroutine genperco

  subroutine stepperco(perco,sca,p)
    implicit none
    integer :: perco(sca),sca
    real :: p,pr
    integer :: i,j,k
    call sr1and()
    call radm(pr)
    do i=1,sca
       if(perco(i)==0) then
          call radm(pr)
          if(pr<p) then
             perco(i)=1
          else
             perco(i)=0
          end if
       end if

    end do
  end subroutine stepperco

