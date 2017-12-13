program integral
implicit none

real(8) :: x, f, x1, f1, s
integer :: i, n

!-----------------------------
s=0.0d0
open(100,file='possbl.dat')
    read(100,'(2x,i8)') n
    read(100,*) x, f
    do i=1,n-1
        read(100,*) x1, f1
        s=s+(f1+f)/2*(x1-x)
        x=x1; f=f1
    enddo
close(100)
!-----------------------------
write(*,*) s


end program integral