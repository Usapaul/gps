program proizvodnaya
implicit none

real(8) :: x, f, x1, f1
integer :: i, n
!-----------------------------
open(100,file='integral.dat')
open(777,file='possbl.dat')
    read(100,'(2x,i8)') n
    write(777,'("# ",i8)') n
    read(100,*) x, f
    do i=1,n
        read(100,*) x1, f1
        write(777,*) x, (f1-f)/(x1-x)
        x=x1; f=f1
    enddo
close(777)
close(100)

end program proizvodnaya